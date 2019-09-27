# in parenthesis
string <- c("(((+X)))", "((+X)-(+Y))")
string %>% 
  gsub("\\(+([^\\(].*?[^\\)])\\)+", "\\1", .)

# string <- c("((<+X>))", "(<+X>-<+Y>)", "(<+X>+<+Y>)/2-<+Z>", "<+X>/2+<+Y>/2-<+Z>")
string <- c("<-KO.1> -    WT", "(<-KO.1> + <+KI.1>)/2 - WT", "(<-KO.1 + KI.1>) - (<+KI.1> + <-KO.1> + <+KI.1>)/3")
string %>%
  gsub("^.*\\[(.*)\\].*", "\\1", .) %>% # may have random terms at the end
  gsub("\\\"", "", .) %>% 
  str_split(",\\s*", simplify = TRUE) %>% 
  gsub("\\s+", "", .) %>% 
  gsub("<([^>]*?)\\+([^>]*?)>", "<\\1.plus.\\2>", .) %>% 
  gsub("<([^>]*?)\\-([^>]*?)>", "<\\1.minus.\\2>", .) %>% 
  gsub("[ <>]+", "", .)


