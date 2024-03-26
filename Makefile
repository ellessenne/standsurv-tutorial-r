.PHONY: style everything

style:
	R -e "styler::style_dir(filetype = c('r', 'rmd', 'qmd'))"

everything:
	make style
	R -f "01-data.R"
	R -f "02-km.R"
	R -f "03-hr.R"
	R -f "04-surv-adj.R"
	R -f "05-surv-std.R"
