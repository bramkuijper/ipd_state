library("tidyverse")
library("patchwork")

the_data <- read_delim(file="summary.csv",delim=";") %>%
    mutate(
        game_type_label=recode(as_factor(game_type),"0"="PD","1"="Snowdrift")
        ) %>%
    rename(Invest_plasticity=xp,Choice_plasticity=yp,Invest=x,Choice=y)


trait_values_l <- the_data %>% select(startup_cost, game_type_label, Invest, Choice) %>%
    pivot_longer(
    cols=c(Invest,Choice)
    ,names_to="trait"
    ,values_to="trait_values")
    
plasticity_values_l <- the_data %>% 
    select(startup_cost, game_type_label, Invest_plasticity, Choice_plasticity) %>%
    pivot_longer(
        cols = c(Invest_plasticity,Choice_plasticity)
        ,names_to="plasticity"
        ,values_to="plasticity_values"
    )

resource_values_l <- the_data %>% 
    select(startup_cost, game_type_label, mean_resources, mean_resources_paired, mean_resources_single) %>%
    pivot_longer(
        cols = c(mean_resources, mean_resources_paired, mean_resources_single)
        ,names_to="resources"
        ,values_to="resource_values"
    )

p0 <- ggplot(data=trait_values_l
       ,mapping=aes(x=startup_cost, y=trait_values)) +
    geom_point(mapping=aes(colour=trait)) +
    facet_grid(.~game_type_label) +
    theme_classic() +
    scale_colour_brewer(palette = "Set1")

p1 <-ggplot(data=plasticity_values_l
       ,mapping=aes(x=startup_cost, y=plasticity_values)) +
    geom_point(mapping=aes(colour=plasticity)) +
    facet_grid(.~game_type_label) +
    theme_classic() +
    scale_colour_brewer(palette = "Set1")

p2 <-ggplot(data=resource_values_l
       ,mapping=aes(x=startup_cost, y=resource_values)) +
    geom_point(mapping=aes(colour=resource_values_ls)) +
    facet_grid(.~game_type_label) +
    theme_classic() +
    scale_colour_brewer(palette = "Set1")

p0/p1/p2

ggsave(filename = "overview_plot.pdf")
