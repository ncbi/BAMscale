#!/usr/bin/env bash

export SCPATH=$(dirname "$BASH_SOURCE")
Rscript -e 'library(methods); shiny::runApp(paste0(Sys.getenv("SCPATH"), "/plotXY/PlotXY.R"), launch.browser = TRUE)'
