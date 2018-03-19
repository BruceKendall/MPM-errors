library(ProjectTemplate)
load.project()

if (!exists("first.stage.comadre")) {
  first_state_comadre <- transmute(flat_comadre,
                                   classnames_first = get_element(classnames))
  write.csv(sort(unique(first_state_comadre$classnames_first)), 
            file = "data/first_stage_comadre.csv", 
            quote = FALSE, row.names = FALSE)
  stop(print("Go edit data/first_stage_comadre.csv"))
}

pre_states <- filter(first.state.comadre, classification == "Pre") %>%
  select(classnames_first) %>% unlist()
post_states <- filter(first.state.comadre, classification == "Post") %>%
  select(classnames_first) %>% unlist()
flat_comadre <- mutate(flat_comadre,
                       first_stage = get_element(classnames),
                       census = ifelse(first_stage %in% post_states, "post", 
                                       ifelse(first_stage %in% pre_states, "pre", "unk"))
                       )

