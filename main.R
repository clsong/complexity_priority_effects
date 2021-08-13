source("toolbox.r")

# sample assembly graphs -------------------------------------------------
total_species <- 3
nodes <- generate_node(LETTERS[1:total_species]) %>%
  group_by(num_species) %>%
  mutate(
    y = num_species,
    x = n()
  ) %>%
  mutate(x = if_else(x > 1, seq(from = 0, to = max(x), length.out = x), total_species / 2)) %>%
  ungroup()

plan(multisession) #parallel
all_possibilities <- nodes %>%
  rename(from = node) %>%
  mutate(invasion = future_map(from, ~ setdiff(LETTERS[1:total_species], str_split(., ",")[[1]]))) %>%
  unnest(invasion) %>%
  select(-num_species) %>%
  mutate(to = future_map2(from, invasion, ~ generate_to(.x, .y))) %>%
  unnest(to) %>%
  mutate(links = paste(from, to, invasion, sep = "-")) %>%
  group_split(from, invasion) %>%
  map(~ pull(., links))

Nsample <- 500
edges <- 1:Nsample %>%
  future_map_dfr(function(x) {
    sample_list <- all_possibilities %>%
      map(~ sampleWithoutSurprises(.))
    names(sample_list) <- paste0("Var", 1:length(sample_list))
    sample_list %>%
      as_tibble()
  }) %>% 
  mutate(graph_label = row_number()) %>%
  gather(key, value, -graph_label) %>%
  group_split(graph_label) %>%
  future_map(~ separate(., value, into = c("from", "to", "invasion"), sep = "-"), .progress = T) %>%
  future_map(~ select(., from, to, invasion))

graphs <- edges %>%
  future_map(~ tbl_graph(nodes, .), .progress = T)

# plot network ------------------------------------------------------------
plot_transition_graph(graphs[[10]])


# classification ----------------------------------------------------------
stable_states <- graphs %>%
  future_map(~ determine_stable_states(.), .progress = T)

transient_path <- stable_states %>%
  map_dbl(~ length(.)) %>%
  enframe(
    name = "graph_label",
    value = "stable_state_num"
  ) %>%
  mutate(stable_states = stable_states) %>%
  mutate(graphs = graphs) %>%
  group_by(stable_state_num) %>%
  mutate(transient_path = future_map2(graphs, stable_states, ~ determine_transient_path(.x, .y), .progress = T)) %>%
  mutate(transient_path_num = future_map(transient_path, function(x) {
    if (nrow(x) > 0) {
      x %>%
        distinct() %>%
        count(from, to) %>%
        pull(n) %>% 
        max()
    } else {
      0
    }
  }, .progress =T)) %>%
  unnest(transient_path_num)

composition_cycle <- transient_path %>% 
  mutate(cycles = future_map(graphs, ~FindCycles(.), .progress = T)) %>% 
  mutate(cycles_num = future_map(cycles, ~length(.))) %>% 
  unnest(cycles_num)

df <- composition_cycle %>%
  ungroup() %>%
  mutate(overlap_between_stable_cycles = future_pmap(list(stable_states, cycles, graphs), 
                                                     determine_overlap_between_stable_cycles)) %>%
  unnest(overlap_between_stable_cycles) %>%
  mutate(overlap_between_stable_cycles = as.factor(overlap_between_stable_cycles)) %>% 
  mutate(cycle_length = map(cycles, function(x) max(unlist(map(x, ~length(.)))))) %>% 
  unnest(cycle_length) %>% 
  select("graph_label", "graphs", "stable_states", "transient_path", "cycles", "stable_state_num", 
         "transient_path_num",  "cycles_num", "cycle_length", "overlap_between_stable_cycles"
  )

# predictability ----------------------------------------------------------
added_edge <- nodes %>%
  select(from = node) %>%
  add_row(from = paste(LETTERS[1:total_species], collapse = ",")) %>%
  unique() %>%
  mutate(invasion = future_map(from, ~ str_split(., ",")[[1]])) %>%
  unnest(invasion) %>%
  mutate(to = from)

ntimes_all <- c(1,5)
entropy2 <- list()
for(ntimes in ntimes_all){
  print(ntimes)
  total_sequence <- rep(LETTERS[1:total_species], ntimes)
  given_sequence <- 1:50 %>%
    map(~sample(total_sequence))
  entropy_random <- edges %>%
    future_map(~generate_community_sequence(., invasion_type = "given", given_sequence = given_sequence), .progress = T) %>%
    future_map(~generate_end_state(., invasion_type = "once"), .progress = T) %>%
    future_map_dbl(~calculate_entropy(.), .progress = T)
  entropy2[[ntimes]] <- entropy_random %>%
    enframe(name = "graph_label", value = "entropy") %>%
    mutate(ntimes = ntimes)
}

df_w_entropy <- df %>% 
  left_join(
    entropy2 %>% 
      bind_rows() %>% 
      nest(entropy = c(entropy, ntimes)),
    by = "graph_label"
  ) %>% 
  unnest(entropy) %>% 
  mutate(cycle_length = ifelse(is.finite(cycle_length), cycle_length, 0))

# prediction --------------------------------------------------------------
data <- df_w_entropy %>%
  select(-graph_label, -stable_states, -transient_path, -cycles, -graphs) %>%
  mutate_all(as.numeric)

library(keras)

fit_neural_network <- function(data, ntime){
  data_longer <- data %>% 
    filter(ntimes == ntime)
  
  # set.seed(1010)
  train_index <- sample(1:nrow(data_longer), 0.8 * nrow(data_longer))
  test_index <- setdiff(1:nrow(data_longer), train_index)
  X <- data_longer %>% 
    mutate(cycles_num = if_else(overlap_between_stable_cycles == 0, 0, cycles_num)) %>%
    select(-entropy)
  y <- data_longer %>% 
    pull(entropy)
  
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  
  model <- keras_model_sequential() %>%
    layer_dense(
      units = 64, activation = "relu",
      input_shape = dim(X_train)[2]
    ) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = 1) %>%
    compile(
      loss = "mse",
      optimizer = optimizer_rmsprop(),
      metrics = list("mean_absolute_error")
    )
  
  model %>% 
    fit(
      x = as.matrix(X_train), y = y_train,
      epochs = 30,
      validation_split = 0.3,
      verbose = 1
    )
  
  predict(model, as.matrix(X_test)) %>% 
    as_tibble() %>% 
    rename(predicted = V1) %>% 
    mutate(observed = y_test) %>% 
    mutate(ntimes = ntime)
}

prediction <- ntimes_all %>% 
  map(~fit_neural_network(data, .)) %>% 
  bind_rows()

prediction_predictability <- prediction %>% 
  mutate(
    predicted = 1- predicted,
    observed = 1-observed
  ) 

# figure -----------------------------------------------------------------
p1 <- prediction_predictability %>%
  filter(ntimes == 1) %>%
  ggplot(aes(predicted, observed)) +
  ggpointdensity::geom_pointdensity() +
  viridis::scale_color_viridis() +
  geom_abline(slope = 1, intercept = 0, color = "gray67", size = .8) +
  jtools::theme_nice() +
  labs(
    x = "Observed predictability",
    y = "Classified predictability",
    title = "Invasion time = 1"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  coord_fixed() +
  theme(
    text = element_text(size = 20),
    legend.position = c(.85, .17),
    axis.line = element_line(color = "black", size = 1),
    legend.title = element_blank()
  )

p2 <- prediction %>%
  nest(values = -ntimes) %>%
  mutate(cor = map(values, ~ cor.test(.$predicted, .$observed))) %>% 
  mutate(
    estimate = map(cor, ~.$estimate),
    confidence_lower = map(cor, ~.$conf.int[1]),
    confidence_upper = map(cor, ~.$conf.int[2])
  ) %>% 
  unnest(estimate, confidence_lower, confidence_upper) %>% 
  ggplot(aes(ntimes, estimate)) +
  scale_x_continuous(breaks = ntimes_all) +
  ylim(0.5, 1) +
  geom_point(size = 4) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = confidence_lower, ymax = confidence_upper), width = .1) +
  labs(
    x = "Invasion times",
    y = "Correlation between\nobserved and classified",
    title = "Explanatory power"
  ) +
  jtools::theme_nice() +
  theme(
    text = element_text(size = 20),
    legend.position = c(.85, .17),
    axis.line = element_line(color = "black", size = 1),
    legend.title = element_blank()
  )

p3 <-
  prediction_predictability %>%
  filter(ntimes == max(ntimes_all)) %>%
  ggplot(aes(predicted, observed)) +
  ggpointdensity::geom_pointdensity() +
  viridis::scale_color_viridis() +
  geom_abline(slope = 1, intercept = 0, color = "gray67", size = .8) +
  jtools::theme_nice() +
  labs(
    x = "Observed predictability",
    y = "Classified predictability",
    title = "Invasion times = 10"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  coord_fixed() +
  theme(
    text = element_text(size = 20),
    legend.position = c(.85, .17),
    axis.line = element_line(color = "black", size = 1),
    legend.title = element_blank()
  )
library("broom")
library(directlabels)

feature_importance <-  ntimes_all %>%
  map(function(x) {
    lm(entropy ~ stable_state_num + transient_path_num + cycle_length +
         factor(overlap_between_stable_cycles), data = filter(data, ntimes == x))
  }) %>%
  map(function(x) {
    relaimpo::calc.relimp(x)$lmg %>%
      enframe()
  }) %>%
  bind_rows(.id = "ntimes") %>%
  nest(data = -ntimes) %>% 
  mutate(ntimes = ntimes_all) %>%
  unnest(data) %>% 
  mutate(
    name = case_when(
      name == "stable_state_num" ~ "# stable\nstates",
      name == "transient_path_num" ~ "# transient\npaths",
      name == "cycle_length" ~ "Length of\ncycle",
      name == "factor(overlap_between_stable_cycles)" ~ "Exit from cycles\nto stable states"
    )
  ) %>%
  mutate(
    name = ordered(name,
                   levels = c("# stable\nstates", "# transient\npaths", "Length of\ncycle", "Exit from cycles\nto stable states")
    )
  ) %>%
  group_by(ntimes) %>%
  mutate(value = value / sum(value)) %>%
  ungroup()
# mutate(label = if_else(ntimes == 8, as.character(name), NA_character_))

# xposition <- 15 #5species
xposition <- 13
p4 <-
  feature_importance %>% 
  ggplot(aes(ntimes, value, fill = name, group = name)) +
  geom_area(alpha = 0.6, size = 1, colour = "black") +
  scale_fill_manual(values = c("#662D91", "#FF7BAC", "#39B54A", "gray20")) +
  # geom_text(aes(label = name), position = position_stack(vjust = 0.5)) +
  labs(
    x = "Invasion times",
    y = "Relative importance\nto classification",
    title = "Feature importance"
  ) +
  # ggrepel::geom_label_repel(aes(label = label),
  #                           nudge_x = 1,
  #                           na.rm = TRUE
  # ) +
  annotate("label", x = xposition, y = .1, label = "Exit from cycles\nto stable states", size = 4, fill = 'gray', alpha = .2) +
  annotate("label", x = xposition, y = .45, label = "Length of cycle", size = 4, fill = '#39B54A', alpha = .6) +
  annotate("label", x = xposition, y = .62, label = "# transient paths", size = 4,  fill = '#FF7BAC', alpha = .6) +
  annotate("label", x = xposition, y = .955, label = "# stable states", size = 4, fill = '#662D91', alpha = .6, color = 'white') +
  # coord_cartesian(clip = "off") +
  # geom_dl(aes(label = name), method = list(dl.combine("last.points")), cex = 0.9)+
  scale_x_continuous(breaks = ntimes_all) +
  jtools::theme_nice() +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none"
  )

library(patchwork)
p1 + p3 + p2 + p4 + 
  plot_annotation(tag_levels = "A") & theme(plot.title = element_text(size = 22), aspect.ratio = 1)


