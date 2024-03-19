.PHONY: style

style:
	R -e "styler::style_dir(filetype = c('r', 'rmd', 'qmd'))"
