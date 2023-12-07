format_signif <- function(x, digits) {
  formatC(signif(x, digits = digits), digits = digits, format = "fg", flag = "#")
}

format_pval <- function(p, digits) {
  ifelse(p < 10^-digits,
         sprintf("<%s", format(10^-digits, digits = 1, scientific = FALSE)),
         sprintf(paste0("%0.", digits, "f"), p))
}
