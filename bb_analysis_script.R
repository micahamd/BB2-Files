# Load packages # -------------------
library('here')
library('tidyverse')
library('TOSTER')
library('rstatix')
library('effectsize')
library('lme4')
library('ARTool')
library('corrplot')
library('pwr')
library('ggthemes')
library('ggpubr')

#-----------------------
# File paths # ---------
study1.path <- "C:/Users/.../bb_study1_data.csv" # Set file path to Experiment 1 CSV file
study2.path <- "C:/Users/.../bb_study2_data.csv" # Set file path to Experiment 2 CSV file

#-----------------------
### Study 1 ### -------
#-----------------------
# Import and clean raw data -----------

bb1.wide <- read.csv(file = study1.path,header = TRUE, na.strings = "NaN") %>% 
  mutate(
    ID = as.factor(studentId),
    Group = as.factor(mcq_answer_1),
    Sex = as.factor(mcq_answer_2),
    Age = as.numeric(sa_answer_11),
    Hearing_Pref = as.factor(mcq_answer_3),
    pre.happy   = as.numeric(sa_answer_1),
    pre.stress  = as.numeric(sa_answer_2),
    pre.peace   = as.numeric(sa_answer_3),
    pre.focus   = as.numeric(sa_answer_4),
    pre.content = as.numeric(sa_answer_5),
    post.happy   = as.numeric(sa_answer_6),
    post.stress  = as.numeric(sa_answer_7),
    post.peace   = as.numeric(sa_answer_8),
    post.focus   = as.numeric(sa_answer_9),
    post.content = as.numeric(sa_answer_10)
  ) # Keep separate for running correlations

bb1.df <- bb1.wide%>%
  mutate(      # Compute normalized change scores
    Peacefulness  = ((post.peace-pre.peace)/(post.peace+pre.peace)), 
    Calmness      = ((post.stress-pre.stress)/(post.stress+pre.stress))*-1, # Reverse score
    Happiness     = ((post.happy-pre.happy)/(post.happy+pre.happy))*-1, # Reverse score
    Focus         = ((post.focus-pre.focus)/(post.focus+pre.focus)), 
    Contentment   = ((post.content-pre.content)/(post.content+pre.content)) 
  ) %>% 
  select(c(ID, Age, Sex, Group, Peacefulness, Calmness, Happiness, Focus, Contentment)) %>% 
  drop_na()%>% 
  mutate(Group = fct_recode(Group,
                            "9 Hz (Alpha-1)" = "NINE HERTZ",
                            "6 Hz (Theta)" = "SIX HERTZ",
                            "3 Hz (Delta)" = "THREE HERTZ",
                            "12 Hz (Alpha-2)" = "TWELVE HERTZ")) %>%
  mutate(Group = fct_relevel(Group, "3 Hz (Delta)", "6 Hz (Theta)", "9 Hz (Alpha-1)", "12 Hz (Alpha-2)"))

# Participant characteristics (S1) ----
bb1.participants <- bb1.df %>%
  group_by(Sex) %>%
  summarize(
    count = n(),
    meanAge = mean(Age, na.rm = TRUE),
    sdAge = sd(Age, na.rm = TRUE),
   .groups = "drop") 

# Convert wide to long data frame for factorial analysis ---------
bb1.long <- bb1.df %>% 
  pivot_longer(
    cols = c(Peacefulness, Calmness, Happiness, Focus, Contentment),
    names_to = "Mood",
    values_to = "change.scores"
  ) %>%
  mutate(Mood = as.factor(Mood))


# Correlations across dependent measures (exploratory)-----------
bb1.cor <- bb1.wide %>%
  select(41:50) %>%
  cor() 
  
# Visualize Correlation Matrix
corrplot(bb1.cor, 
         method = "number",
         diag = FALSE,insig = 'blank',
         col = c("red", "brown", "orange", "yellow"), type = "upper")

# Calculate coefficients range
bb1.coefs <- bb1.cor %>%
  as.data.frame() %>%
  rownames_to_column("Mood1") %>%
  pivot_longer(-Mood1, names_to = "Mood2", values_to = "Coefficients") %>%
  filter(Coefficients != 1)

coef.range <- range(bb1.coefs$Coefficients)








# Create model object using the 'lme4' package  ----
aov.mod <- lme4::lmer(change.scores ~ Group * Mood + (1 | ID), data = bb1.long)

# Model assumption tests ----------
aov.res <- resid(aov.mod)           # Extract Model Residuals
aov.shap <- shapiro.test(aov.res)   # Shapiro Tests on Residuals
aov.cook <- cooks.distance(aov.mod) # Cook's distance for identifying outliers
cook.outliers <- which(aov.cook > 4 / length(aov.res)) # Identify outliers

# Create data frame *without outliers*  ---------
bb1.long.no.outliers <- bb1.long[-cook.outliers, ] # Create data frame without outliers

# Homogeneity of variance -----------
bb1.lev.p.value <- round((rstatix::levene_test(change.scores ~ Group * Mood, data = bb1.long))$p, 3)
bb1.lev.p.value.no.outliers <- round((rstatix::levene_test(change.scores ~ Group * Mood, data = bb1.long.no.outliers))$p, 3)

#-----------------------
#---------- Helper Functions ----------
# For reporting TOST 'run.tost()': ----
run.tost <-function(data){
  df <- data %>% 
    group_by(Group, Mood) %>%
    get_summary_stats(change.scores,type = "full") %>%
    mutate(low.ci = mean - (ci / 2), up.ci = mean + (ci / 2)) %>%
    mutate(
      equiv.bound = case_when(
        # 1. Test for rejection (CI is completely outside the bounds)
        low.ci > eqb_value | up.ci < -eqb_value ~ "REJECTED",
        
        # 2. Test for acceptance (CI is completely inside the bounds)
        low.ci > -eqb_value & up.ci < eqb_value ~ "ACCEPTED",
        
        # 3. All other cases are undecided (the CI overlaps one or both bounds)
        TRUE                            ~ "UNDECIDED"
      ),
      ROPE = as.factor(equiv.bound)
    ) %>%
    rowwise() %>%
    # Run tsum_TOST() for each row.
    mutate(tost_output = list(
      tsum_TOST(m1 = mean, n1 = n, sd1 = sd, eqb = eqb_value, mu = 0)
    )) %>%
    mutate(
      # Hedges's g values
      g_estimate = tost_output$effsize["Hedges's g", "estimate"],
      g_lower_ci = tost_output$effsize["Hedges's g", "lower.ci"],
      g_upper_ci = tost_output$effsize["Hedges's g", "upper.ci"],
      
      # EQ bounds
      raw_eqb_lo = tost_output$eqb$low_eq[1],
      raw_eqb_up = tost_output$eqb$high_eq[1],
      g_eqb_lo = tost_output$eqb$low_eq[2],
      g_eqb_up = tost_output$eqb$high_eq[2],
      
      # Decision criteria
      eq_test = tost_output$decision$TOST,
      ot_test = tost_output$decision$ttest
    ) %>%
    ungroup() %>%
    select(-c(tost_output, variable, ci, min, max, q1, q3, iqr, mad,raw_eqb_lo,raw_eqb_up,equiv.bound))
  df
} # Run and output TOSTs with bias corrected effects
# For running Type 3 ANOVA on ART transformed data 'run_art_anova(data)': ----
run_art_anova <- function(data, suffix = "") {
  art.mod <- ARTool::art(change.scores ~ Group * Mood + (1|ID), data = data)
  aov_res <- anova(art.mod, response = "art") %>% as.data.frame()
  art.effects <- as.data.frame((effectsize::cohens_f(art.mod))[,1:2])
  aov_res <- full_join(aov_res,art.effects,by = c("Term"="Parameter"))
  
  col_name <- paste0("res.art", suffix)
  
  aov_res %>%
    as.data.frame() %>%
    mutate(!!col_name := paste0(
      "F-ART", if(suffix != "") "/NoOutliers" else "", " (", Df, ", ", round(Df.res), ") = ", round(F, 3),
      ", p = ", ifelse(`Pr(>F)` > .001, round(`Pr(>F)`, 3), "0.001"), ", f = ", round(Cohens_f_partial,2)
    )) %>%
    mutate(Term = as.factor(Term)) %>%
    select(Term, !!col_name)
}     # Type 3 ANOVA after ART
# For running Type 3 ANOVA normally 'run_rstatix_anova(data)': ----
run_rstatix_anova <- function(data, suffix = "") {
  aov_res <- rstatix::anova_test(
    data = data, between = Group, within = Mood,
    wid = ID, dv = change.scores, type = 3, effect.size = "pes",
    white.adjust = TRUE
  )
  col_name <- paste0("res.no.art", suffix)
  
  aov_res$ANOVA %>%
    as.data.frame() %>%
    mutate(!!col_name := paste0(
      "F-NoART", if(suffix != "") "/NoOutliers" else "", " (", DFn, ", ", round(DFd), ") = ", round(F, 3),
      ", p = ", ifelse(p > .001, round(p, 3), "0.001"),
      ", Cohen's f = ", round(effectsize::eta2_to_f(pes), 3)
    )) %>%
    mutate(Term = as.factor(Effect)) %>%
    select(Term, !!col_name)
} # Type 3 ANOVA on rstatix
# For drawing TOST plots 'tost.plot(data,xlabel="",ylabel="")': -----
tost.plot <- function(data,xlabel="",ylabel=""){
  
  xlab_text <- (if (xlabel == "") "Group" else xlabel)
  ylab_text <- (if (ylabel == "") "Mean change scores with 95% CIs" else ylabel)
  
  data %>% 
    ggplot(aes(x = Group, y = mean, color = ROPE)) +
    geom_rect(
      mapping = aes(xmin = -Inf, xmax = Inf, ymin = -eqb_value, ymax = eqb_value),
      fill = "lightgray", color = "white", alpha = 0.1
    )+
    geom_pointrange(aes(ymin = low.ci, ymax = up.ci), size = .4, alpha = .9, position = position_dodge(.6)) +
    labs(x = xlab_text,y =  ylab_text, fill = "Mood")+
    ylim(-.2,.4)+
    geom_hline(yintercept = 0, linetype = "dotted") +
    facet_wrap(~ Mood, nrow = 1)+
    ggthemes::theme_few() +
    theme(legend.position = "top")
}

# For drawing effect size plots 'g.plot(data,xlabel="",ylabel="")': -----
g.plot <- function(data,xlabel="",ylabel=""){
  
  xlab_text <- (if (xlabel == "") "Group" else xlabel)
  ylab_text <- (if (ylabel == "") "Hedge's g scores with 95% CIs" else ylabel)
  
  data %>% 
  ggplot(aes(x = Group, y = g_estimate, color = ROPE)) +
    geom_pointrange(aes(ymin = g_lower_ci, ymax = g_upper_ci), size = .8, alpha = .9, position = position_dodge(.6)) +
    labs(
      x = xlab_text, y = ylab_text, fill = "Mood"
    ) +
    ylim(-.5,1.5)+
    geom_hline(yintercept = 0, linetype = "dotted") +
    facet_wrap(~ Mood, nrow = 1) +
    ggthemes::theme_few() +
    theme(legend.position = "none")
  
}
# For drawing descriptives plot 'ind.plot(data,ylabel="")': -----------
ind.plot<-function(data,ylabel=""){
    ggplot(data,aes(x=Mood, y=change.scores, color = Mood,fill=Mood))+
    geom_violin(alpha=.3,color=NA)+
    # geom_boxplot(alpha=.3)+
    geom_point(position = position_dodge2(width=.4),size=1.4,alpha=.9)+
    facet_wrap(~Group)+
    ylab(ylabel)+
    geom_hline(yintercept=0,linetype="dashed")+
    theme(legend.position="none")+
    ggthemes::theme_few()
}  

#-----------------------
# TOST Equivalence Testing # ------
#-----------------------
# Find smallest effect that can be determined with 95% power for the smallest sample -------
smallest.sample <- bb1.df %>% group_by(Group) %>% summarise(sample = n()) %>% pull() %>% min()
smallest.effect <- pwr::pwr.t.test(n=smallest.sample, sig.level=0.05, power=0.8, type="one.sample", alternative = "greater")$d

# Define Equivalence bound -----------

# Identify pooled SD across the sample and multiply with smallest effect of interest
eqb_value <- sd(bb1.long$change.scores)*smallest.effect

# Descriptive and TOST summaries for outlier retained and outlier removed data -----------
bb1.tost <- run.tost(bb1.long)
bb1.no.outliers.tost <- run.tost(bb1.long.no.outliers)

# Plot change score summaries for outlier retained and outlier removed data: 'p1.all.tost.plots' -----------
p1.tost.plot <- tost.plot(bb1.tost,ylabel = " ",xlabel="Study 1\n(Full dataset)") + coord_flip()
p1.no.outliers.tost.plot <- tost.plot(bb1.no.outliers.tost,ylabel = " ",xlabel="Study 1\n(Outliers removed)") + coord_flip()
p1.all.tost.plots <- ggpubr::ggarrange(p1.tost.plot,p1.no.outliers.tost.plot, nrow = 2, labels = c("A.", "B."), common.legend = T)

# Plot Hedge's g estimates for outlier retained and outlier removed data -----------
g1.plot <- g.plot(bb1.tost,ylabel = "Hedge's g with\n95% CIs (full dataset)")
g1.no.outliers.plot <- g.plot(bb1.no.outliers.tost,ylabel = "Hedge's g with\n95% CIs (outliers removed)")
g1.all.plots <- ggpubr::ggarrange(g1.plot,g1.no.outliers.plot, ncol = 2, labels = c("C.", "D."), common.legend = T)

# Combine Plots? 
# ggpubr::ggarrange(p1.all.tost.plots,g1.all.plots,nrow=2)

#-----------------------
# Analyses of Variance (S1) # --------
#-----------------------
# Type 3 ART-ANOVAs # --------
art_full <- run_art_anova(bb1.long, suffix = ".with.out")
art_no_outliers <- run_art_anova(bb1.long.no.outliers, suffix = ".no.out")

# Type 3 ANOVA as-is --------
aov_full <- run_rstatix_anova(bb1.long)
aov_no_outliers <- run_rstatix_anova(bb1.long.no.outliers, suffix = ".no.out")

# Combine ANOVA results into a single data frame 'bb1.all.aov' -----------
bb1.all.aov <- list(art_full, art_no_outliers, aov_full, aov_no_outliers) %>% reduce(full_join, by = "Term")
# Post-hoc Tukey tests ## -------
bb1.tukey <- aov(formula = change.scores ~ Group * Mood, data = bb1.long) %>% tukey_hsd() %>% filter(p.adj.signif != "ns")  # On Full dataset
bb1.no.outlier.tukey <- aov(formula = change.scores ~ Group * Mood, data = bb1.long.no.outliers) %>% tukey_hsd() %>% filter(p.adj.signif != "ns")  # On dataset with outliers removed
#-----------------------
## Output tables for Study 1 ## -----

# bb1.tost %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 1/BB Study 1 Output Tables/bb1.full.descriptives.csv")
# bb1.no.outliers.tost %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 1/BB Study 1 Output Tables/bb1.no.outliers.descriptives.csv")
# bb1.tukey %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 1/BB Study 1 Output Tables/bb1.full.tukey.csv")
# bb1.no.outlier.tukey %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 1/BB Study 1 Output Tables/bb1.no.outliers.tukey.csv")
# bb1.all.aov %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 1/BB Study 1 Output Tables/bb1.all.aov.csv")


#-----------------------
### Study 2 ### -------
#-----------------------
# Import data and compute normalized scores  --------
bb2.df <- read.csv(file = study2.path,header = TRUE, na.strings = "NaN") %>% 
  as.data.frame() %>% 
  mutate(ID = as.factor(studentId),
       pre.happy   = as.numeric(sa_answer_1),
       pre.stress  = as.numeric(sa_answer_2),
       pre.peace   = as.numeric(sa_answer_3),
       pre.focus   = as.numeric(sa_answer_4),
       pre.content = as.numeric(sa_answer_5),
       Group = as.factor(mcq_answer_1),  
       post.happy   = as.numeric(sa_answer_6),
       post.stress  = as.numeric(sa_answer_7),
       post.peace   = as.numeric(sa_answer_8),
       post.focus   = as.numeric(sa_answer_9),
       post.content = as.numeric(sa_answer_10),
       Age = as.numeric(sa_answer_11),
       Sex = as.factor(mcq_answer_2),
       subjective_experience = as.character(Q1_answer),
       comfort_level = as.character(Q2_answer),
       willing_to_repeat = as.character(Q3_answer),
       subjective_experience.coded = as.factor(Q1_rate),
       comfort_level.coded = as.factor(Q2_rate),
       willing_to_repeat.coded = as.factor(Q3_rate)
       ) %>% 
  mutate(      # Compute normalized change scores
    Peacefulness   = ((post.peace-pre.peace)/(post.peace+pre.peace)), 
      Calmness     = ((post.stress-pre.stress)/(post.stress+pre.stress))*-1,
      Happiness    = ((post.happy-pre.happy)/(post.happy+pre.happy)), 
      Focus        = ((post.focus-pre.focus)/(post.focus+pre.focus)), 
      Contentment  = ((post.content-pre.content)/(post.content+pre.content)) 
  ) %>% 
  select(c(ID, Age, Sex, Group, Peacefulness, Calmness, Happiness, Focus, Contentment, 
           subjective_experience, comfort_level, willing_to_repeat,
           subjective_experience.coded, comfort_level.coded, willing_to_repeat.coded))

# Participant characteristics (S2) ----
bb2.participants <- bb2.df %>%
  group_by(Group) %>%
  summarize(
    count = n(),
    meanAge = mean(Age, na.rm = TRUE),
    sdAge = sd(Age, na.rm = TRUE),
    .groups = "drop") 


# Select data for factorial analyses and convert to long ----------
bb2.long <- bb2.df %>% 
  select(c(ID, Age, Group, Peacefulness, Calmness, Happiness, Focus, Contentment)) %>% 
  pivot_longer(
    cols = c(Peacefulness, Calmness, Happiness, Focus, Contentment),
    names_to = "Mood",
    values_to = "change.scores"
  ) %>%
  mutate(Mood = as.factor(Mood))  

# Create model object using the 'lme4' package  ----
aov.mod2 <- lme4::lmer(change.scores ~ Group * Mood + (1 | ID), data = bb2.long)

# Model assumption tests ----------
aov.res2 <- resid(aov.mod2)           # Extract Model Residuals
aov.shap2 <- shapiro.test(aov.res2)   # Shapiro Tests on Residuals
aov.cook2 <- cooks.distance(aov.mod2) # Cook's distance for identifying outliers
cook.outliers2 <- which(aov.cook2 > 4 / length(aov.res2)) # Identify outliers

# Histograms, QQ plots and Cook's distance visualization ---------
par(mfrow = c(2, 3))

hist(aov.res, main = "Histogram of Residuals", xlab = "Residuals",ylab="Study 1")
qqnorm(aov.res)
qqline(aov.res)

plot(aov.cook, main = "Cook's Distance Plot") # Plot Cook's distance
abline(h = 4 / length(residuals(aov.mod)), col = "red")

hist(aov.res2, main = "Histogram of Residuals", xlab = "Residuals",ylab="Study 2")
qqnorm(aov.res2)
qqline(aov.res2)

plot(aov.cook, main = "Cook's Distance Plot") # Plot Cook's distance
abline(h = 4 / length(residuals(aov.mod)), col = "red")


# Create data frame *without outliers*  ---------
bb2.long.no.outliers <- bb2.long[-cook.outliers2, ] # Create data frame without outliers


# Homogeneity of variance -----------
bb2.lev.p.value <- round((rstatix::levene_test(change.scores ~ Group * Mood, data = bb2.long))$p, 3)
bb2.lev.p.value.no.outliers <- round((rstatix::levene_test(change.scores ~ Group * Mood, data = bb2.long.no.outliers))$p, 3)

# Descriptive and TOST summaries for outlier retained and outlier removed data -----------
bb2.tost <- run.tost(bb2.long)
bb2.no.outliers.tost <- run.tost(bb2.long.no.outliers)


#-----------------------
# Analyses of Variance (S2) # --------
#-----------------------
# Type 3 ART-ANOVAs # --------
art_full <- run_art_anova(bb2.long, suffix = ".with.out")
art_no_outliers <- run_art_anova(bb2.long.no.outliers, suffix = ".no.out")

# Type 3 ANOVA as-is --------
aov_full <- run_rstatix_anova(bb2.long)
aov_no_outliers <- run_rstatix_anova(bb2.long.no.outliers, suffix = ".no.out")

# Combine ANOVA results into a single data frame 'bb2.all.aov' -----------
bb2.all.aov <- list(art_full, art_no_outliers, aov_full, aov_no_outliers) %>% reduce(full_join, by = "Term")

# Post-hoc Tukey tests ## -------
bb2.tukey <- aov(formula = change.scores ~ Group * Mood, data = bb2.long) %>% tukey_hsd() %>% filter(p.adj.signif != "ns")  # On Full dataset
bb2.no.outlier.tukey <- aov(formula = change.scores ~ Group * Mood, data = bb2.long.no.outliers) %>% tukey_hsd() %>% filter(p.adj.signif != "ns")  # On dataset with outliers removed


# Plot effects for study 2 ---------
p2.tost.plot <- tost.plot(bb2.tost,ylabel = " ",xlabel="Study 2\n(Full dataset)") + coord_flip()
p2.no.outliers.tost.plot <- tost.plot(bb2.no.outliers.tost,ylabel = "Mean change scores with 95% CIs",xlabel="Study 2\n(Outliers removed)") + coord_flip()
p2.all.tost.plots <- ggpubr::ggarrange(p2.tost.plot,p2.no.outliers.tost.plot, nrow = 2, labels = c("C.", "D."), common.legend = T)

# Combine all plots from both studies ---------
# ggpubr::ggarrange(p1.all.tost.plots,p2.all.tost.plots, nrow=2,common.legend = T) # Everything
ggpubr::ggarrange(p1.no.outliers.tost.plot,p2.no.outliers.tost.plot, nrow=2,common.legend = T) # Outliers removed dataset

## Output tables for Study 2 ## -----

# bb2.tost %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 2/bb2.full.descriptives.csv")
# bb2.no.outliers.tost %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 2/bb2.no.outliers.descriptives.csv")
# bb2.tukey %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 2/bb2.full.tukey.csv")
# bb2.no.outlier.tukey %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 2/bb2.no.outliers.tukey.csv")
# bb2.all.aov %>% write.csv(file = "C:/Users/micah/OneDrive/Desktop/Papers/binaural-begin-v2/binaural manuscript/Exp 2/bb2.all.aov.csv")

#-----------------------
# Sentiment Analyses ---------
#-----------------------
# Correlation across sentiments (coded manually) ---------
qual.cor <- bb2.df %>% 
  mutate(
  "Subjective Experience"= as.factor(subjective_experience.coded),
  "Comfort Level"= as.factor(comfort_level.coded),
  "Willing to Repeat"= as.factor(willing_to_repeat.coded)
  ) %>% 
  select(`Subjective Experience`,`Comfort Level`,`Willing to Repeat`,
         Peacefulness,Calmness, Happiness, Focus,Contentment) %>% 
  mutate_all(as.numeric) %>% 
  scale() %>% 
  cor(method = "spearman") 

# Calculate p-values from correlation matrix
qual.p.mat <- (corrplot::cor.mtest(qual.cor))$p %>% as.matrix()
qual.cor.plot <- corrplot::corrplot(qual.cor,
                   type = "upper",
                   diag = F,
                   sig.level=.05,
                   p.mat = qual.p.mat,
                   method = 'number',
                   pch.cex = 1,
                   pch.col="grey",
                   insig = "pch")


# (corrplot::cor.mtest(qual.cor))$p # p.values

# Examine frequency distributions of coded responses -----------------
bb2.freq <-bb2.df %>% 
  select(Group,Sex,subjective_experience.coded,comfort_level.coded,willing_to_repeat.coded) %>% 
  pivot_longer(cols=c(subjective_experience.coded,comfort_level.coded,willing_to_repeat.coded),names_to = "Qual.Category",values_to = "Coded.Responses") %>% 
  mutate(Qual.Category = as.factor(Qual.Category))%>%
  mutate(Qual.Category=fct_recode(Qual.Category,
                                  "Subjective Experience" = "subjective_experience.coded",
                                  "Comfort Level" = "comfort_level.coded",
                                  "Willing to Repeat" = "willing_to_repeat.coded")) %>% 
  mutate(Coded.Responses=fct_recode(Coded.Responses,
                                    "Positive"="1",
                                    "Uncertain"="0",
                                    "Negative"="-1")) %>% 
  mutate(`Qualitative Category` = as.factor(Qual.Category)) %>% 
  select(-Qual.Category) 

# Explore frequency percentages ----

bb2.freq %>%
  count(Group,`Qualitative Category`,Coded.Responses) %>% 
  group_by(Group,`Qualitative Category`) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup() %>% 
  filter(Coded.Responses=="Negative")
  


# Plot frequency distributions of coded responses -----------------
bb2.freq.plot <- bb2.freq %>% 
  ggplot(aes(x = Coded.Responses, fill = `Qualitative Category`)) +
  geom_bar(position = position_dodge(width = 0.6), width = 0.5, alpha = 0.9) +
  geom_text(stat = "count",
            aes(label = ..count..),
            position = position_dodge(width = 0.6),
            vjust = -0.3, size = 4) +
  ylab("Frequency") + xlab("Coded Sentiments") +
  ylim(0,30)+
  facet_wrap(~ Group) +
  ggthemes::theme_base()+
  scale_fill_brewer(palette = "Blues")+
  theme(legend.position="top")


# Chi-square tests

# Create the contingency table
bb2.freq.tab <- bb2.freq %>%
  count(Group, Coded.Responses) %>%
  pivot_wider(names_from = Coded.Responses, values_from = n, values_fill = 0) %>%
  column_to_rownames("Group") %>%
  as.matrix()

# Run the chi-square test
chisq.test(bb2.freq.tab[1,])
chisq.test(bb2.freq.tab[2,])
chisq.test(bb2.freq.tab[3,])
chisq.test(bb2.freq.tab[4,])

# ---------------------
# Descriptive plots -------------

i1.p <- ind.plot(bb1.long,ylabel="Normalized Change Scores (Study 1)")
i2.p <- ind.plot(bb2.long,ylabel="Normalized Change Scores (Study 2)")
i.all.p <- ggpubr::ggarrange(i1.p,i2.p,labels = c("A","B"),nrow=2,common.legend = T)




# ---------------------