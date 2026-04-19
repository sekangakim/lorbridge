## data-raw/create_datasets.R
## Run once to generate the .rda files in data/
## usethis::use_data() saves each object to data/<name>.rda

vm_raw_txt <- "
54 0 0 1 0
59 0 1 0 0
62 0 1 0 0
63 4 0 2 0
65 0 0 1 0
67 2 2 1 0
69 2 4 1 0
71 2 1 0 0
73 1 2 1 0
74 4 3 2 1
76 7 5 1 1
78 6 3 2 1
80 4 5 2 0
81 9 4 2 0
82 7 2 6 1
84 10 4 6 0
85 9 6 7 0
86 4 5 7 0
87 18 8 10 2
89 13 6 2 2
90 14 6 3 1
92 18 4 5 2
93 23 5 3 2
95 16 7 6 1
96 22 10 2 3
98 34 6 4 4
100 27 3 8 0
101 22 2 3 2
103 22 6 4 4
104 16 3 4 2
105 16 6 7 4
107 29 5 8 5
108 22 2 5 1
110 21 0 5 0
112 18 2 3 1
113 14 3 1 4
115 19 1 2 0
117 14 1 5 0
119 11 0 3 2
121 20 1 2 0
123 13 0 2 3
125 9 0 2 1
126 9 1 1 1
128 9 0 1 0
130 7 0 1 1
132 3 2 0 1
134 4 0 0 0
136 2 0 1 0
138 2 0 0 0
139 1 1 1 0
143 0 1 0 0
147 0 0 1 0
149 1 0 0 0
"

vm_raw <- read.table(text = vm_raw_txt, header = FALSE,
                     col.names = c("VM", "Race4", "Race1", "Race2", "Race3"))

vm_breaks <- c(-Inf, 64.28, 81.00, 100.00, 121.00, 138.36, Inf)
vm_labels <- paste0("VM", 1:6)
vm_raw$VMbin <- cut(vm_raw$VM, breaks = vm_breaks,
                    labels = vm_labels, right = TRUE, include.lowest = TRUE)

# Expand to individual-level data
lorbridge_data <- do.call(rbind, lapply(seq_len(nrow(vm_raw)), function(i) {
  vm_val  <- vm_raw$VM[i]
  bin_val <- as.character(vm_raw$VMbin[i])
  race_cols     <- c("Race4", "Race1", "Race2", "Race3")
  minority_flag <- c(Race4 = 0L, Race1 = 1L, Race2 = 1L, Race3 = 1L)
  do.call(rbind, lapply(race_cols, function(r) {
    n <- vm_raw[[r]][i]
    if (n == 0L) return(NULL)
    data.frame(VM       = rep(vm_val, n),
               VMbin    = rep(bin_val, n),
               minority = rep(minority_flag[r], n),
               Race     = rep(r, n),
               stringsAsFactors = FALSE)
  }))
}))
row.names(lorbridge_data) <- NULL
lorbridge_data$VMbin <- relevel(factor(lorbridge_data$VMbin, levels = vm_labels),
                                 ref = "VM4")

# IQ contingency table (4 races x 6 IQ bins)
tab_IQ <- matrix(
  c(30, 42, 32, 25,  9,  2,
    28, 38, 38, 24, 10,  9,
     2,  7, 10, 12, 11, 11,
    15, 51,121,172,125, 76),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("Race1","Race2","Race3","Race4"), paste0("IQ", 1:6))
)

# VM contingency table (4 races x 6 VM bins)
tab_VM <- matrix(
  c( 2, 29,  72,  32, 3, 2,
     3, 13,  69,  52, 8, 2,
     0,  3,  18,  25, 7, 0,
     4, 37, 215, 244,58, 2),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("Race1","Race2","Race3","Race4"), paste0("VM", 1:6))
)

# IQ x VM table (6 x 6, for DONSCA)
tab_IQ_VM <- matrix(
  c(28,30,14, 3, 0, 0,
    15,59,36,20, 7, 1,
    10,44,66,47,23,11,
    10,20,59,88,33,23,
     1, 7,19,55,43,30,
     1, 1, 7,25,24,40),
  nrow = 6, byrow = TRUE,
  dimnames = list(paste0("IQ", 1:6), paste0("VM", 1:6))
)

usethis::use_data(lorbridge_data, tab_IQ, tab_VM, tab_IQ_VM,
                  overwrite = TRUE)
