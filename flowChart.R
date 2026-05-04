# initialize connection to DB ---------------------------------------------
library(DBI)
dbcon <- dbConnect(RPostgres::Postgres(),
                   user = "gavriel.fialkoff@senseerahealth.com",
                   password = system(paste0("gcloud sql generate-login-token"), 
                                     intern = TRUE),
                   host = "34.165.66.84",
                   dbname = "senseera-prod",
                   sslmode='require')
dbGetQuery(dbcon, "SELECT *
                    FROM sample_qc_view;") %>% 
  mutate(sample_name_short = sub("_.*", "", sample_name)) -> all.qc
sample.annotation %>% 
  count(group)

all.qc %>% 
  filter(grepl("SCLC", sample_name)) %>% 
  group_by(sample_name_short) %>% 
  slice_max(order_by = tss_pct_signal) %>% 
  filter(tss_pct_signal < 25 | tss_win_count < 5e5) %>%
  filter(!sample_name_short %in% sample.annotation$Sample_id) %>% 
  select(sample_name_short, tss_pct_signal, tss_win_count) %>% 
  dim()
  

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(DiagrammeR)

p <- grViz("
digraph flow_chart {

  # Global settings for nodes
  node [fontname = Helvetica, 
        shape = rectangle, 
        style = filled, 
        fillcolor = White, 
        penwidth = 1.5, 
        fontsize = 7] # Smaller font size

  # Define Nodes (removed bold style)
  node1 [label = 'SCLC/NEC plasma samples\n(n=353)']
  node2 [label = 'cfChIP-seq samples\n(n=285)']
  node3 [label = 'control groups (Healthy/CRC/NSCLC)\n(n=156)', fillcolor = '#fafafa']
  node4 [label = 'final cohort\n(n=441)', fillcolor = '#f0f8ff']

  # Define Edges
  edge [penwidth = 1, arrowsize = 0.8]
  node1 -> node2 [label = ' quality filtering', fontsize = 7]
  node2 -> node4
  node3 -> node4

  {rank = same; node2; node3;}
}
")
p
# Exporting with high resolution
p %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg_png(
    file = paste0(figDirPaper, "Study_Flow_Chart.png"),
    width = 2400,  # Increased for higher density
    height = 1800  # Adjust height to match your desired aspect ratio
  )
