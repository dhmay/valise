# methods for calculating sensitivity and specificity, and building ROC curves

calc_sens_spec <- function(predictions, outcome, higher_is_better=T) {
  isCanc = outcome[order(predictions,decreasing=higher_is_better)]
  sensitivity <- cumsum(isCanc)/sum(isCanc)
  specificity <- cumsum(1-isCanc)/sum(1-isCanc)  
  cbind(sensitivity, specificity)
}

calc_auc <- function(sensitivity, specificity) { sum(diff(specificity)*(sensitivity[-1]+sensitivity[-length(sensitivity)]))/2 }

# nicely-formatted ROC curve
plot_roc <- function(sensitivity, specificity, color, newplot, title, linetype) {
	if (newplot) {
		plot(specificity, sensitivity, type="l", xlab="1-Specificity",ylab="Sensitivity", col=color, lwd=2,cex=2, main=title, lty=linetype)
	}
	else {
		lines(specificity, sensitivity, col=color, lwd=2, lty=linetype)
	}
	abline(0,1,col = 'gray',lty = 2)
	abline(v=axis(1),col = 'gray',lty = 2)
	abline(h=axis(2),col = 'gray',lty = 2)
}

#formatting for legend rows
buildLegendRow<-function(label, auc) {
  paste(label,'. AUC=',round(auc,digits=2), sep='');
}

calc_plot_roc <- function(predictions, outcome, color, newplot, title='',linetype=1,
                          show_auc=FALSE,
                          higher_is_better=T) {
    sensAndSpec = calc_sens_spec(predictions, outcome, higher_is_better=higher_is_better)
    auc = calc_auc(sensAndSpec[,1], sensAndSpec[,2])
    if (show_auc) {
      title = paste(title,' (AUC=',round(auc,3),')', sep='')
    }
    plot_roc(sensAndSpec[,1], sensAndSpec[,2], color, newplot, title,linetype) 
    calc_auc(sensAndSpec[,1], sensAndSpec[,2])
    auc
}

# interpolate
find_sens_at_spec <- function(specs, senss, spec) {
  ind = 0
  for (i in 1:length(specs)) {
    if (specs[i] < spec) {
      ind = i
      break
    }
  }
  proportion = (spec - specs[ind-1]) / (specs[ind] - specs[ind-1])
  senss[ind-1] + proportion * (senss[ind] - senss[ind-1])
}

