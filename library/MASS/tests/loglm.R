###  Tests to exercise more of loglm()
## Added for 7.3-60
library(MASS)

minn38_tab <- xtabs(f ~ ., minn38)

mod_llm <- loglm(~ hs*phs*fol + sex, minn38_tab)      ## log-linear
mod_glm <- glm(f ~ hs*phs*fol + sex, poisson, minn38) ## surrogate Poisson

## deviance check
all.equal(deviance(mod_llm), deviance(mod_glm))

## residuals check
minn38_res <- within(minn38, {
  res_glm_d <- resid(mod_glm)
  res_glm_p <- resid(mod_glm, type = "pearson")
  res_glm_r <- resid(mod_glm, type = "response")
})

.tmp <- capture.output({
  .tmp_d <- as.table(resid(mod_llm))                          ## array to table
  .tmp_d <- as.data.frame(.tmp_d, responseName = "res_llm_d") ## table to data frame

  .tmp_p <- as.table(resid(mod_llm, type = "pearson"))        ## array to table
  .tmp_p <- as.data.frame(.tmp_p, responseName = "res_llm_p") ## table to data frame

  .tmp_r <- as.table(resid(mod_llm, type = "response"))       ## array to table
  .tmp_r <- as.data.frame(.tmp_r, responseName = "res_llm_r") ## table to data frame
})

minn38_res <- merge(merge(merge(minn38_res, .tmp_d), .tmp_p), .tmp_r)

with(minn38_res, c(
     "deviance OK" = isTRUE(all.equal(res_llm_d, res_glm_d)),
     "pearson OK"  = isTRUE(all.equal(res_llm_p, res_glm_p)),
     "resoibse OK" = isTRUE(all.equal(res_llm_r, res_glm_r))))

cat(.tmp, sep = "\n")                         ## need to allow this to be quiet

## update and anova check

mod_llm_2 <- update(mod_llm, . ~ . - sex + sex*hs)
mod_glm_2 <- update(mod_glm, . ~ . - sex + sex*hs)

anova(mod_llm, mod_llm_2, test = "LR")     ### check visually ("LR" could not be "LRT")
anova(mod_glm, mod_glm_2, test = "LR")     ### check visually ("LR" could     be "LRT")

## extractAIC check

daic_llm <- extractAIC(mod_llm_2) - extractAIC(mod_llm)
daic_glm <- extractAIC(mod_glm_2) - extractAIC(mod_glm)

all.equal(daic_llm, daic_glm)

## Alternative check
d_llm <- dropterm(mod_llm_2, test = "Chisq")
d_glm <- dropterm(mod_glm_2, test = "Chisq")

all.equal(diff(d_llm$AIC), diff(d_glm$AIC))

# rm(.tmp, .tmp_d, .tmp_p, .tmp_r)
# rm(list = c("d_glm", "d_llm", "daic_glm", "daic_llm", "minn38_res",
#             "minn38_tab", "mod_glm", "mod_glm_2", "mod_llm", "mod_llm_2"))

