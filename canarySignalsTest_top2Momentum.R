

# ============================================================
# Canary Signal Test
# What should tell a tactical strategy to reduce risk?
# ============================================================

library(quantmod)
library(PerformanceAnalytics)
library(xts)
library(zoo)

stratStats <- function(rets, digits = 4) {
  
  rets <- na.omit(rets)
  annual <- apply.yearly(rets, Return.cumulative)
  
  stats <- rbind(
    Return.annualized(rets),
    StdDev.annualized(rets),
    SharpeRatio.annualized(rets, Rf = 0),
    maxDrawdown(rets),
    CalmarRatio(rets),
    apply(annual, 2, min, na.rm = TRUE)
  )
  
  rownames(stats) <- c(
    "Annualized Return",
    "Annualized Volatility",
    "Annualized Sharpe",
    "Worst Drawdown",
    "Calmar Ratio",
    "Worst Calendar Year"
  )
  
  return(round(stats, digits))
}




run_canary_signal_test <- function(dataStartDate = "2006-01-01",
                                   analysisStartDate = "2007-07-01",
                                   endDate = Sys.Date(),
                                   investableAssets = c("SPY", "EFA", "EEM", "VNQ", "DBC"),
                                   canaryPairs = list(
                                     VWO_BND = c("VWO", "BND"),
                                     EFA_BND = c("EFA", "BND"),
                                     EEM_BND = c("EEM", "BND"),
                                     HYG_BND = c("HYG", "BND"),
                                     DBC_BND = c("DBC", "BND"),
                                     SPY_BND = c("SPY", "BND")
                                   ),
                                   nAssets = 2,
                                   assetMomWindows = c(63, 126),
                                   canaryMomWindows = c(20, 60),
                                   riskOffMode = c("conditional_defensive", "defensive", "cash"),
                                   defensiveAsset = "IEF",
                                   defensiveMomWindows = c(20, 60),
                                   cashAsset = NULL,
                                   cashReturn = 0,
                                   rebalanceOn = "months",
                                   verbose = TRUE) {
  
  riskOffMode <- match.arg(riskOffMode)
  
  # -----------------------------
  # Helper functions
  # -----------------------------
  cumulative_return <- function(x) {
    prod(1 + as.numeric(x), na.rm = FALSE) - 1
  }
  
  momentum_score <- function(R, symbol, endRow, windows) {
    
    scores <- sapply(windows, function(w) {
      
      startRow <- endRow - w + 1
      
      if (startRow < 1) {
        return(NA_real_)
      }
      
      cumulative_return(R[startRow:endRow, symbol])
    })
    
    sum(scores, na.rm = FALSE)
  }
  
  # -----------------------------
  # Symbols needed
  # -----------------------------
  canarySymbols <- unique(unlist(canaryPairs))
  
  downloadSymbols <- unique(c(
    investableAssets,
    canarySymbols,
    defensiveAsset,
    cashAsset
  ))
  
  downloadSymbols <- downloadSymbols[!is.na(downloadSymbols)]
  
  if (is.null(cashAsset)) {
    cashName <- "CASH"
  } else {
    cashName <- cashAsset
  }
  
  tradeAssets <- unique(c(
    investableAssets,
    defensiveAsset,
    cashName
  ))
  
  # -----------------------------
  # Download adjusted prices
  # -----------------------------
  priceList <- list()
  
  for (sym in downloadSymbols) {
    
    if (verbose) {
      message("Downloading: ", sym)
    }
    
    tmp <- getSymbols(
      sym,
      from = dataStartDate,
      to = endDate,
      auto.assign = FALSE,
      warnings = FALSE
    )
    
    px <- Ad(tmp)
    colnames(px) <- sym
    priceList[[sym]] <- px
  }
  
  prices <- do.call(merge, priceList)
  prices <- na.omit(prices)
  
  assetReturns <- Return.calculate(prices)
  assetReturns <- na.omit(assetReturns)
  
  if (is.null(cashAsset)) {
    
    cashReturns <- xts(
      rep(cashReturn / 252, NROW(assetReturns)),
      order.by = index(assetReturns)
    )
    
    colnames(cashReturns) <- "CASH"
    
    allReturns <- merge(assetReturns, cashReturns)
    
  } else {
    
    allReturns <- assetReturns
  }
  
  allReturns <- na.omit(allReturns)
  
  if (verbose) {
    message("Return data begins: ", as.character(first(index(allReturns))))
    message("Return data ends:   ", as.character(last(index(allReturns))))
  }
  
  # -----------------------------
  # Rebalance dates
  # -----------------------------
  rebalanceRows <- endpoints(allReturns, on = rebalanceOn)
  rebalanceRows <- rebalanceRows[rebalanceRows > 0]
  rebalanceRows <- rebalanceRows[rebalanceRows <= NROW(allReturns)]
  
  maxLookback <- max(
    assetMomWindows,
    canaryMomWindows,
    defensiveMomWindows
  )
  
  signalRows <- rebalanceRows[rebalanceRows > maxLookback]
  
  # ============================================================
  # Baseline 1: Equal Weight
  # ============================================================
  equalWeight <- Return.portfolio(
    R = allReturns[, investableAssets],
    weights = rep(1 / length(investableAssets), length(investableAssets)),
    rebalance_on = rebalanceOn
  )
  
  colnames(equalWeight) <- "EqualWeight"
  
  # ============================================================
  # Baseline 2: Top-N Momentum without Canary
  # ============================================================
  run_topn_momentum <- function() {
    
    weightsMat <- matrix(
      NA_real_,
      nrow = NROW(allReturns),
      ncol = length(investableAssets)
    )
    
    colnames(weightsMat) <- investableAssets
    
    for (signalRow in signalRows) {
      
      assetScores <- sapply(investableAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, assetMomWindows)
      })
      
      rankedAssets <- names(sort(assetScores, decreasing = TRUE))
      selectedAssets <- rankedAssets[1:nAssets]
      
      w <- rep(0, length(investableAssets))
      names(w) <- investableAssets
      w[selectedAssets] <- 1 / nAssets
      
      weightsMat[signalRow, ] <- w[investableAssets]
    }
    
    weights <- xts(weightsMat, order.by = index(allReturns))
    weights <- na.locf(weights, na.rm = FALSE)
    weightsLag <- lag(weights, k = 1)
    
    rets <- xts(
      rowSums(weightsLag * allReturns[, investableAssets], na.rm = FALSE),
      order.by = index(allReturns)
    )
    
    colnames(rets) <- paste0("Top", nAssets, "_Momentum")
    
    return(list(
      returns = rets,
      weights = weights,
      weightsLag = weightsLag
    ))
  }
  
  topN <- run_topn_momentum()
  topNReturns <- topN$returns
  
  # ============================================================
  # Canary strategy engine
  # ============================================================
  run_one_canary <- function(pairName, canaryAssets) {
    
    weightsMat <- matrix(
      NA_real_,
      nrow = NROW(allReturns),
      ncol = length(tradeAssets)
    )
    
    colnames(weightsMat) <- tradeAssets
    
    riskBudgetVec <- rep(NA_real_, NROW(allReturns))
    
    signalRecords <- vector("list", length(signalRows))
    recordCounter <- 1
    
    for (signalRow in signalRows) {
      
      signalDate <- index(allReturns)[signalRow]
      
      # -----------------------------
      # Rank investable assets by momentum
      # -----------------------------
      assetScores <- sapply(investableAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, assetMomWindows)
      })
      
      rankedAssets <- names(sort(assetScores, decreasing = TRUE))
      selectedAssets <- rankedAssets[1:nAssets]
      
      # -----------------------------
      # Canary risk budget
      # -----------------------------
      canaryScores <- sapply(canaryAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, canaryMomWindows)
      })
      
      riskBudget <- mean(canaryScores > 0)
      riskOffWeight <- 1 - riskBudget
      
      # -----------------------------
      # Start weights
      # -----------------------------
      w <- rep(0, length(tradeAssets))
      names(w) <- tradeAssets
      
      # Risk-on allocation
      w[selectedAssets] <- riskBudget / nAssets
      
      # Risk-off allocation
      if (riskOffWeight > 0) {
        
        if (riskOffMode == "cash") {
          
          w[cashName] <- w[cashName] + riskOffWeight
          
        } else if (riskOffMode == "defensive") {
          
          w[defensiveAsset] <- w[defensiveAsset] + riskOffWeight
          
        } else if (riskOffMode == "conditional_defensive") {
          
          defensiveScore <- momentum_score(
            allReturns,
            defensiveAsset,
            signalRow,
            defensiveMomWindows
          )
          
          if (!is.na(defensiveScore) && defensiveScore > 0) {
            w[defensiveAsset] <- w[defensiveAsset] + riskOffWeight
          } else {
            w[cashName] <- w[cashName] + riskOffWeight
          }
        }
      }
      
      weightsMat[signalRow, ] <- w[tradeAssets]
      riskBudgetVec[signalRow] <- riskBudget
      
      signalRecords[[recordCounter]] <- data.frame(
        date = as.Date(signalDate),
        canaryPair = pairName,
        canaryAssets = paste(canaryAssets, collapse = ", "),
        selectedAssets = paste(selectedAssets, collapse = ", "),
        riskBudget = riskBudget,
        riskOffWeight = riskOffWeight,
        stringsAsFactors = FALSE
      )
      
      recordCounter <- recordCounter + 1
    }
    
    # -----------------------------
    # Convert signal weights to daily weights
    # -----------------------------
    weights <- xts(weightsMat, order.by = index(allReturns))
    weights <- na.locf(weights, na.rm = FALSE)
    weightsLag <- lag(weights, k = 1)
    
    # Strategy returns
    rets <- xts(
      rowSums(weightsLag * allReturns[, tradeAssets], na.rm = FALSE),
      order.by = index(allReturns)
    )
    
    colnames(rets) <- pairName
    
    # Live period
    retsLive <- rets[paste0(analysisStartDate, "/")]
    retsLive <- na.omit(retsLive)
    
    liveWeights <- weightsLag[index(retsLive), tradeAssets]
    
    riskBudget <- xts(riskBudgetVec, order.by = index(allReturns))
    riskBudget <- na.locf(riskBudget, na.rm = FALSE)
    riskBudgetLag <- lag(riskBudget, k = 1)
    liveRiskBudget <- riskBudgetLag[index(retsLive)]
    
    signalLog <- do.call(rbind, signalRecords)
    
    # -----------------------------
    # Exposure stats
    # -----------------------------
    exposure <- data.frame(
      Strategy = pairName,
      avgRiskBudget = mean(as.numeric(liveRiskBudget), na.rm = TRUE),
      pctFullRiskOn = mean(as.numeric(liveRiskBudget) == 1, na.rm = TRUE),
      pctHalfRiskOn = mean(as.numeric(liveRiskBudget) == 0.5, na.rm = TRUE),
      pctRiskOff = mean(as.numeric(liveRiskBudget) == 0, na.rm = TRUE),
      avgRiskAssetWeight = mean(rowSums(liveWeights[, investableAssets], na.rm = TRUE), na.rm = TRUE),
      avgDefensiveWeight = mean(as.numeric(liveWeights[, defensiveAsset]), na.rm = TRUE),
      avgCashWeight = mean(as.numeric(liveWeights[, cashName]), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    exposure[, -1] <- round(exposure[, -1], 4)
    
    return(list(
      returns = retsLive,
      weights = weights,
      weightsLag = weightsLag,
      riskBudget = riskBudget,
      signalLog = signalLog,
      exposure = exposure
    ))
  }
  
  # -----------------------------
  # Run all canary pairs
  # -----------------------------
  canaryResults <- list()
  canaryReturnsList <- list()
  exposureList <- list()
  
  for (pairName in names(canaryPairs)) {
    
    if (verbose) {
      message("Running canary pair: ", pairName)
    }
    
    res <- run_one_canary(
      pairName = pairName,
      canaryAssets = canaryPairs[[pairName]]
    )
    
    canaryResults[[pairName]] <- res
    canaryReturnsList[[pairName]] <- res$returns
    exposureList[[pairName]] <- res$exposure
  }
  
  canaryReturns <- do.call(merge, canaryReturnsList)
  exposureStats <- do.call(rbind, exposureList)
  rownames(exposureStats) <- NULL
  
  # -----------------------------
  # Combine all returns
  # -----------------------------
  mainReturns <- merge(
    equalWeight,
    topNReturns,
    canaryReturns
  )
  
  mainReturns <- mainReturns[paste0(analysisStartDate, "/")]
  mainReturns <- na.omit(mainReturns)
  
  canaryReturns <- canaryReturns[index(mainReturns)]
  canaryReturns <- na.omit(canaryReturns)
  
  # -----------------------------
  # Summary tables
  # -----------------------------
  summary <- stratStats(mainReturns)
  canarySummary <- stratStats(canaryReturns)
  
  summaryTable <- data.frame(
    Strategy = colnames(mainReturns),
    AnnualizedReturn = as.numeric(summary["Annualized Return", ]),
    AnnualizedVolatility = as.numeric(summary["Annualized Volatility", ]),
    Sharpe = as.numeric(summary["Annualized Sharpe", ]),
    WorstDrawdown = as.numeric(summary["Worst Drawdown", ]),
    Calmar = as.numeric(summary["Calmar Ratio", ]),
    WorstCalendarYear = as.numeric(summary["Worst Calendar Year", ]),
    stringsAsFactors = FALSE
  )
  
  rankedBySharpe <- summaryTable[order(-summaryTable$Sharpe), ]
  rankedByCalmar <- summaryTable[order(-summaryTable$Calmar), ]
  rankedByDrawdown <- summaryTable[order(summaryTable$WorstDrawdown), ]
  
  rownames(rankedBySharpe) <- NULL
  rownames(rankedByCalmar) <- NULL
  rownames(rankedByDrawdown) <- NULL
  
  annualReturns <- apply.yearly(mainReturns, Return.cumulative)
  
  return(list(
    mainReturns = mainReturns,
    canaryReturns = canaryReturns,
    summary = summary,
    canarySummary = canarySummary,
    summaryTable = summaryTable,
    rankedBySharpe = rankedBySharpe,
    rankedByCalmar = rankedByCalmar,
    rankedByDrawdown = rankedByDrawdown,
    annualReturns = round(annualReturns, 4),
    exposureStats = exposureStats,
    canaryResults = canaryResults,
    prices = prices,
    allReturns = allReturns,
    settings = list(
      dataStartDate = dataStartDate,
      analysisStartDate = as.character(first(index(mainReturns))),
      endDate = as.character(last(index(mainReturns))),
      investableAssets = investableAssets,
      canaryPairs = canaryPairs,
      nAssets = nAssets,
      assetMomWindows = assetMomWindows,
      canaryMomWindows = canaryMomWindows,
      riskOffMode = riskOffMode,
      defensiveAsset = defensiveAsset,
      defensiveMomWindows = defensiveMomWindows,
      cashAsset = cashAsset,
      cashReturn = cashReturn,
      rebalanceOn = rebalanceOn
    )
  ))
}




canaryTest <- run_canary_signal_test(
  dataStartDate = "2006-01-01",
  analysisStartDate = "2007-07-01",
  endDate = "2026-07-09",
  investableAssets = c("SPY", "EFA", "EEM", "VNQ", "DBC"),
  canaryPairs = list(
    VWO_BND = c("VWO", "BND"),
    EFA_BND = c("EFA", "BND"),
    EEM_BND = c("EEM", "BND"),
    HYG_BND = c("HYG", "BND"),
    DBC_BND = c("DBC", "BND"),
    SPY_BND = c("SPY", "BND")
  ),
  nAssets = 2,
  assetMomWindows = c(63, 126),
  canaryMomWindows = c(20, 60),
  riskOffMode = "conditional_defensive",
  defensiveAsset = "IEF",
  defensiveMomWindows = c(20, 60),
  cashAsset = NULL,
  cashReturn = 0,
  rebalanceOn = "months",
  verbose = TRUE
)




canaryTest$summary
canaryTest$canarySummary
canaryTest$summaryTable
canaryTest$rankedBySharpe
canaryTest$rankedByCalmar
canaryTest$rankedByDrawdown
canaryTest$exposureStats
canaryTest$annualReturns
canaryTest$settings




charts.PerformanceSummary(
  canaryTest$mainReturns,
  main = "Testing Canary Signals in a Top 2 Momentum Strategy",
  wealth.index = T
)




charts.PerformanceSummary(
  canaryTest$canaryReturns,
  main = "Canary Signal Comparison",
  wealth.index = T
)




chart.Drawdown(
  canaryTest$canaryReturns,
  main = "Canary Signal Drawdowns",
  legend.loc = "bottom"
)




table.CalendarReturns(canaryTest$mainReturns)




canaryTest$rankedBySharpe
canaryTest$rankedByCalmar
canaryTest$exposureStats










































































