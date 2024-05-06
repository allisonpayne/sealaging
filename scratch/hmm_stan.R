library(tidyverse)
library(patchwork)

# psi = detectability. Random effect so it can vary by year. 
# psi1 = breeding, psi2 = non-breeding
# psi1 ~ beta(psi1_alpha^2, psi1_beta^2)
# mu_psi1 = (sqrt(10), sqrt(2))
# sd_psi1_alpha = sd_psi1_beta = 0.25
# cor_psi1_alpha_beta = -1
# cor_var_psi1_alpha_beta = sd_psi1_alpha * sd_psi1_beta * cor_psi1_alpha_beta
# sigma_psi1 = (sd_psi1_alpha^2         cor_var_psi1_alpha_beta
#               cor_var_psi1_alpha_beta sd_psi1_alpha^)
# psi1_alpha_beta ~ mvnormal(mu_psi1, sigma_psi1)

# psi1 (Breeding detection) -----------------------------------------------

mu_psi1ab <- sqrt(c(10, 2))
sd_psi1a <- sd_psi1b <- 0.25
cor_psi1ab <- -0.5
cor_var_psi1ab <- cor_psi1ab * sd_psi1a * sd_psi1b
sigma_psi1 <- matrix(c(sd_psi1a^2, cor_var_psi1ab, cor_var_psi1ab, sd_psi1b^2), 
                     nrow = 2)
psi1ab <- mvtnorm::rmvnorm(1e5, mean = mu_psi1ab, sigma = sigma_psi1)
psi1 <- rbeta(1e5, psi1ab[, 1]^2, psi1ab[, 2]^2)

# plot
# Psi1 alpha and beta
psi1ab_plot <- tibble(`psi[1*alpha]` = psi1ab[, 1]^2, 
                      `psi[1*beta]` = psi1ab[, 2]^2) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(value, color = name)) + 
  geom_density(linewidth = 1.5) + 
  scale_color_manual(values = c("cornflowerblue", "firebrick"), 
                     labels = scales::parse_format()) + 
  labs(x = "Parameter", y = "Density") +
  theme_classic() + 
  theme(legend.position = c(0.75, 0.75), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))
# Psi1
psi1_plot <- tibble(psi1 = psi1) %>% 
  ggplot(aes(psi1)) +
  geom_density(linewidth = 1.5) +
  scale_x_continuous(name = expression(psi[1])) +
  labs(y = "Density") +
  theme_classic()

psi1ab_plot / psi1_plot

# psi2 (Non-breeding detection) -------------------------------------------

mu_psi2ab <- sqrt(c(2, 5))
sd_psi2a <- sd_psi2b <- 0.25
cor_psi2ab <- -0.5
cor_var_psi2ab <- cor_psi2ab * sd_psi2a * sd_psi2b
sigma_psi2 <- matrix(c(sd_psi2a^2, cor_var_psi2ab, cor_var_psi2ab, sd_psi2b^2), 
                     nrow = 2)
psi2ab <- mvtnorm::rmvnorm(1e5, mean = mu_psi2ab, sigma = sigma_psi2)
psi2 <- rbeta(1e5, psi2ab[, 1]^2, psi2ab[, 2]^2)

# plot
# Psi1 alpha and beta
psi2ab_plot <- tibble(`psi[2*alpha]` = psi2ab[, 1]^2, 
                      `psi[2*beta]` = psi2ab[, 2]^2) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(value, color = name)) + 
  geom_density(linewidth = 1.5) + 
  scale_color_manual(values = c("cornflowerblue", "firebrick"), 
                     labels = scales::parse_format()) + 
  labs(x = "Parameter", y = "Density") +
  theme_classic() + 
  theme(legend.position = c(0.75, 0.75), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))
# Psi1
psi2_plot <- tibble(psi2 = psi2) %>% 
  ggplot(aes(psi2)) +
  geom_density(linewidth = 1.5) +
  scale_x_continuous(name = expression(psi[2])) +
  labs(y = "Density") +
  theme_classic()

psi2ab_plot / psi2_plot
