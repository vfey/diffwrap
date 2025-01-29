# test package using RALP data

## sample sheet
si.f <- "/Users/fsvife/Workspace/UTU/Dhana/diffr_test/input_example/cleaned_sample_sheet.txt"
samp.info <- read.delim(si.f)
samp.info <- diff_expr_get_samp_info(samp.info, samples = "sample", groups = "group")

## read counts
in.fls <- dir("/Users/fsvife/Workspace/UTU/Dhana/diffr_test/input_example", pattern="^RALP", full=T)
names(in.fls) <- sub("_cleaned", "", basename(in.fls))
counts <- diff_expr_read_counts(in.fls, samp.info = samp.info)

## design matrix
groups <- relevel(samp.info$Groups, ref="normal")
design <- diff_expr_make_design(samp.info, groups, NULL, F)
