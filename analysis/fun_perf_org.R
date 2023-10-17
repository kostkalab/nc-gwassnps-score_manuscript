library(pROC)
library(PRROC)

make_roc_org = function(phenotype, ctrl, case, sc){
  print(phenotype)
  ctrl_p = ctrl[ctrl$phenotype == phenotype]
  case_p = case[case$phenotype == phenotype]
  type = c(ctrl_p$type, case_p$type)
  score = rbind(mcols(ctrl_p)[sc], mcols(case_p)[sc])
  score = as.data.frame(score)

  ## for roc
    auroc = integer(length(score))
    for (i in seq_along(score)){
      auroc[i] = roc(type, score[,i])$auc
    }
    auc_roc = data.frame(PTPE = phenotype, auc = auroc, variant_score = sc)
    auc_roc$measure = "roc"

  ## for pr curve
    aupr = integer(length(score))
    for (i in seq_along(score)){
      s = !is.na(score[,i])
      aupr[i] = pr.curve(scores.class0=score[s,i],weights.class0=type[s])$auc.integral
    }
    auc_pr = data.frame(PTPE = phenotype, auc = aupr, variant_score = sc)
    auc_pr$measure = "pr"
    auc = rbind(auc_roc, auc_pr)
  return(auc)
}

