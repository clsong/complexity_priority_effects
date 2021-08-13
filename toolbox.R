library(here)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(patchwork)
library(furrr)
library(entropy)

generate_composition <- function(all_species_composition, S) {
  combn(all_species_composition, S) %>%
    t() %>%
    as_tibble() %>%
    unite(node, sep = ",") %>%
    pull(node)
}

generate_node <- function(all_species_composition) {
  1:length(all_species_composition) %>%
    map(~ generate_composition(all_species_composition, .)) %>%
    unlist() %>%
    as_tibble() %>%
    rename(node = value) %>%
    mutate(num_species = (str_length(node) + 1) / 2)
}

plot_transition_graph <- function(graph, edge_type = 'arc') {
  if(edge_type == 'arc'){
    graph %>%
      ggraph(layout = "manual", x = x, y = y) +
      geom_edge_arc(
        arrow = arrow(length = unit(5, "mm")),
        end_cap = circle(5, "mm"),
        aes(color = invasion)
      ) +
      geom_edge_loop(
        arrow = arrow(length = unit(5, "mm")),
        end_cap = circle(5, "mm"),
        aes(color = invasion)
      ) +
      geom_node_point(size = 12, colour = "steelblue") +
      geom_node_text(aes(label = node), colour = "white", vjust = 0.4) +
      theme_graph(base_family = "Helvetica") +
      theme(aspect.ratio = 1)
  } else{
    graph %>%
      ggraph(layout = "manual", x = x, y = y) +
      geom_edge_link(
        arrow = arrow(length = unit(5, "mm")),
        end_cap = circle(5, "mm"),
        aes(color = invasion)
      ) +
      geom_edge_loop(
        arrow = arrow(length = unit(5, "mm")),
        end_cap = circle(5, "mm"),
        aes(color = invasion)
      ) +
      geom_node_point(size = 12, colour = "steelblue") +
      geom_node_text(aes(label = node), colour = "white", vjust = 0.4) +
      theme_graph(base_family = "Helvetica") +
      theme(aspect.ratio = 1)
  }
}

shifter <- function(x, n) {
  if (n == 0) {
    x <- x
  } else {
    x <- c(tail(x, -n), head(x, n))
  }
  x
}

generate_to <- function(from, invasion) {
  species <- c(str_split(from, pattern = ",")[[1]], str_split(invasion, pattern = ",")[[1]]) %>%
    sort()

  1:length(species) %>%
    map(~ generate_composition(species, .)) %>%
    unlist()
}

generate_community_sequence <- function(edge, invasion_type, given_sequence = NA) {
  edge <- bind_rows(edge, added_edge)

  generate_given_sequence <- function(edge, sequence) {
    community_sequence <- vector("list", length = length(sequence))
    community_sequence[[1]] <- sequence[1]
    for (i in 1:(length(sequence) - 1)) {
      community_sequence[[i + 1]] <- edge %>%
        filter(from == community_sequence[[i]] & invasion == sequence[i + 1]) %>%
        pull(to)
    }
    community_sequence <- unlist(community_sequence)
    community_sequence
  }

  if (invasion_type == "sequential") {
    sequence <- rep(LETTERS[1:total_species], 10)
    result <- 0:(total_species - 1) %>%
      map(~ generate_given_sequence(edge, shifter(sequence, .)))
  } else if (invasion_type == "once") {
    result <- gtools::permutations(total_species, total_species, LETTERS[1:total_species], repeats = FALSE) %>%
      t() %>%
      split(rep(1:ncol(.), each = nrow(.))) %>%
      map(~ generate_given_sequence(edge, .))
  } else if (invasion_type == "random") {
    result <- sample(LETTERS[1:total_species], 1000, replace=T) %>% 
      matrix(nrow = 100, ncol = 10) %>% 
      split(rep(1:ncol(.), each = nrow(.))) %>% 
      map(~ generate_given_sequence(edge, .))
  } else if (invasion_type == "given") {
    result <- given_sequence %>% 
      map(~ generate_given_sequence(edge, .))
  }
  result
}

generate_community_sequence_092020 <- function(edge, given_sequences) {
  edge <- bind_rows(edge, added_edge)
  
  given_sequences %>% 
    map(function(sequence){
      community_sequence <- vector("list", length = length(sequence))
      community_sequence[[1]] <- sequence[1]
      for (i in 1:(length(sequence) - 1)) {
        community_sequence[[i + 1]] <- edge[edge$from == community_sequence[[i]] & edge$invasion == sequence[i + 1], ]$to
      }
      unlist(community_sequence)
    })
}

calc_pairwise_overlaps <- function(sets) {
  n_sets <- length(sets)
  set_names <- 1:n_sets
  n_overlaps <- choose(n = n_sets, k = 2)
  
  vec_name1 <- character(length = n_overlaps)
  vec_name2 <- character(length = n_overlaps)
  vec_overlap <- integer(length = n_overlaps)
  overlaps_index <- 1
  
  for (i in seq_len(n_sets - 1)) {
    name1 <- set_names[i]
    set1 <- sets[[i]]
    for (j in seq(i + 1, n_sets)) {
      name2 <- set_names[j]
      set2 <- sets[[j]]
      
      vec_name1[overlaps_index] <- name1
      vec_name2[overlaps_index] <- name2
      vec_overlap[overlaps_index] <- if_else(any(all_equal(set1, set2) == TRUE), 1, 0)
      
      overlaps_index <- overlaps_index + 1
    }
  }
  tibble(
    name1 = vec_name1,
    name2 = vec_name2,
    overlap = vec_overlap
  )
}

sampleWithoutSurprises <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}

identify_duplicated_graphs <- function(edges, mutated_rows = c("from", "to", "invasion"), sort = TRUE) {
  permutations <- combinat::permn(1:total_species)
  generate_permuted <- function(edge, mutated_rows, sort) {
    permutations %>%
      map(function(x) {
        edge %>%
          mutate_at(mutated_rows, function(y) textclean::mgsub(y, LETTERS[1:total_species], x)) %>%
          ungroup() %>%
          rowwise() %>%
          mutate_at(mutated_rows, function(y) {
            if(sort == TRUE){
              paste(sort(LETTERS[as.numeric(str_split(y, ",")[[1]])]), collapse = ",")
            } else {
              paste(LETTERS[as.numeric(str_split(y, ",")[[1]])], collapse = ",")
            }
            }) %>%
          ungroup()
      }) %>%
      bind_rows(.id = "permutation_label")
  }

  df <- edges %>%
    map(~ generate_permuted(., mutated_rows, sort)) %>%
    bind_rows(.id = "graph_label") %>%
    nest(data = -c(graph_label, permutation_label))

  duplicated <- calc_pairwise_overlaps(df$data) %>%
    filter(overlap == 1) %>%
    mutate(
      name1 = (as.integer(name1) - 1) %/% factorial(total_species) + 1,
      name2 = (as.integer(name2) - 1) %/% factorial(total_species) + 1
    ) %>%
    unique() %>%
    filter(name1 != name2)

  duplicated
}

unique_graphs <- function(edges_all){
  
  find_characters <- function(edge){
    all_characters <- expand.grid(
      self_loop = 0:1,
      from_num = 1:(total_species-1),
      diff_num = c(-1,0,1) 
    ) 
    all_characters2 <-   expand.grid(
        num_species = 1:total_species,
        m = 0:5
      )
  
    bind_rows(
      edge %>% 
        mutate(
          self_loop = if_else(from == to, 1, 0),
          from_num = (str_length(from) + 1) / 2,
          diff_num = (str_length(to) - str_length(from)) / 2
        ) %>%
        count(., self_loop, from_num, diff_num) %>% 
        full_join(all_characters, by = c("self_loop", "from_num", "diff_num")) %>% 
        arrange(self_loop, from_num, diff_num) %>% 
        select(n)
    ) %>% 
      replace_na(list(n = 0))
  }
  
  graph_groups <- edges_all %>%
    future_map(~find_characters(.), .progress= T) %>% 
    bind_rows(.id = 'graph_label') %>%
    select(graph_label, n) %>%
    group_by(graph_label) %>%
    mutate(characters = row_number()) %>%
    pivot_wider(names_from = characters, values_from = n) %>%
    group_by_at(setdiff(names(.),  "graph_label")) 

  graph_groups_split <- graph_groups %>% 
    group_split() 
  
  identify_duplicated_within_group <- function(graph_group){
    pb$tick()$print()
    graph_group <- graph_group %>% 
      pull(graph_label) %>% 
      as.numeric() %>% 
      enframe() %>% 
      rename(graph_label = value) 
    
    identify_duplicated_graphs(edges_all[graph_group$graph_label]) %>% 
      left_join(graph_group, by = c("name1" = "name")) %>% 
      left_join(graph_group, by = c("name2" = "name")) %>% 
      select(graph_label.x, graph_label.y)
  }
 
  pb <- progress_estimated(length(graph_groups_split))
  
  duplicated_graphs <- graph_groups_split %>% 
    map(~identify_duplicated_within_group(.))
  
  duplicated_graphs
}

unique_invasion_sequence <- function(end_states){
  end_states_summary <- end_states %>%
    select(graph_label, end_state) %>%
    mutate(num_species = (str_length(end_state) + 1) / 2) %>%
    nest(data = -graph_label) %>%
    mutate(data = future_map(data, function(x) {
      bind_rows(
        x %>% count(num_species) %>% rename(total_n = n) %>% gather(key, value, -num_species),
        x %>% distinct() %>% count(num_species) %>% rename(unique_n = n) %>% gather(key, value, -num_species)
      ) %>%
        full_join(.,
                  expand.grid(
                    num_species = 1:total_species,
                    key = c("total_n", "unique_n")
                  ),
                  by = c("num_species", "key")
        ) %>%
        replace_na(list(value = 0))
    })) %>%
    unnest(data) %>%
    mutate(key = paste(key, num_species, sep = "-")) %>%
    select(-num_species) %>%
    pivot_wider(names_from = key, values_from = value) %>%
    group_by_at(setdiff(names(.), "graph_label"))  
  
  groups <- end_states_summary %>% 
    group_split() %>% 
    map(~pull(., graph_label))
  
  df <- end_states %>%
    select(graph_label, end_state) %>%
    group_by(graph_label) %>%
    mutate(invasion_label = row_number()) %>%
    left_join(
      combinat::permn(LETTERS[1:total_species]) %>% 
        map(~paste(., collapse = ',')) %>% 
        unlist() %>% 
        enframe(name = "invasion_label", value = 'invasion_sequence'),
      by = "invasion_label"
    ) %>% 
    select(-invasion_label)
  
  identify_duplicated_within_group <- function(df, group){
    group_name <- enframe(group)
    df %>% 
      filter(graph_label %in% group) %>% 
      group_by(graph_label) %>% 
      group_split() %>% 
      map(~select(., -graph_label)) %>% 
      identify_duplicated_graphs(mutated_rows = c("end_state", "invasion_sequence"), sort = FALSE) %>% 
      left_join(group_name, by = c("name1" = "name")) %>% 
      left_join(group_name, by = c("name2" = "name")) %>% 
      select(value.x, value.y)
  }
  
  groups %>% 
    future_map(~identify_duplicated_within_group(df, .), .progress = T)
}

generate_end_state <- function(pattern, invasion_type) {
  if (invasion_type == "sequential") {
    
  } else if (invasion_type == "once") {
    end_states <- pattern %>%
      map(~ tail(., n = 1)) %>%
      unlist()
    result <- end_states
    # result <- length(unique(end_states))
  }
  result
}

calculate_entropy <- function(x){
  y <- x %>% 
    table() %>% 
    enframe() %>% 
    mutate(value = value/sum(value)) %>% 
    full_join(
      select(nodes, name = node),
      by = "name"
    ) %>% 
    replace_na(list(value = 0)) %>% 
    pull(value)
  
  entropy(y)/log(length(y))
}

determine_overlap_between_stable_cycles <- function(stable_states, cycles, graph){
  stable_node <- stable_states %>% 
    {which(. == nodes$node)}
  
  cycle_node <- cycles %>% 
    unlist() %>% 
    {nodes$node[.]}
  
  if(length(stable_node) > 0 & length(cycle_node) > 0){
    overlap <- graph %>% 
      select(node, num_species) %>%
      mutate(distance = node_distance_to(stable_node)) %>% 
      filter(node %in% cycle_node) %>% 
      filter(!is.infinite(distance)) %>%
      as_tibble() %>%
      nrow()
  } else{
    overlap <- 0
  }
  
  if_else(overlap > 0, 1, 0)
}

determine_stable_states <- function(graph) {
  setdiff(
    1:nrow(nodes) %>%
      map_dbl(function(x) {
        graph %>%
          select(node, num_species) %>%
          mutate(distance = node_distance_from(x)) %>%
          filter(node != nodes$node[x]) %>%
          filter(!is.infinite(distance)) %>%
          as_tibble() %>%
          nrow()
      }) %>% {
        nodes$node[which(. == 0)]
      },
    1:nrow(nodes) %>%
      map_dbl(function(x) {
        graph %>%
          select(node, num_species) %>%
          mutate(distance = node_distance_to(x)) %>%
          filter(num_species == 1) %>%
          filter(!is.infinite(distance)) %>%
          as_tibble() %>%
          nrow()
      }) %>% {
        nodes$node[which(. == 0)]
      }
  )
}

determine_transient_path <- function(graph, stable_state) {
  expand.grid(
    from = 1:total_species,
    to = which(nodes$node %in% stable_state)
  ) %>%
    mutate(transient_path = map2(from, to, ~ igraph::all_simple_paths(graph, from = .x, to = .y, mode = c("out")))) %>%
    mutate(transient_path = as.character(transient_path)) %>%
    as_tibble() %>%
    mutate(transient_path = str_split(transient_path, pattern = ", c")) %>%
    unnest(transient_path) %>%
    mutate(transient_path = map(transient_path, ~ paste(unlist(str_match_all(., "[0-9]+")), collapse = ","))) %>%
    unnest(transient_path)
}

## More efficient version
FindCycles <- function(g) {
  #https://stackoverflow.com/questions/55091438/r-igraph-find-all-cycles
  Cycles = NULL
  num_node <- g %>% as_tibble() %>% nrow()
  # for(v1 in igraph::V(g)) {
  for(v1 in 1:num_node) {
    if(igraph::degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = igraph::neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    connected <- g %>%
      select(node, num_species) %>%
      mutate(distance = node_distance_to(v1)) %>%
      filter(num_species == 1) %>% 
      filter(!is.infinite(distance)) %>%
      as_tibble() %>%
      nrow()
    if(connected > 0){
      for(v2 in GoodNeighbors) {
        TempCyc = lapply(igraph::all_simple_paths(g, v2, v1, mode="out"), function(p) c(v1,p))
        TempCyc = TempCyc[which(sapply(TempCyc, length) > 1)]
        TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
        Cycles  = c(Cycles, TempCyc)
      }
    }
  }
  Cycles
}

