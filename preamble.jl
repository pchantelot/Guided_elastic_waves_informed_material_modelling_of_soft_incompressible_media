# Theme setup
begin
	# Make a plotting theme
	My_theme = Theme(
		fontsize = 16,
		fonts=(; regular=texfont(:text),
                        bold=texfont(:bold),
                        italic=texfont(:italic),
                        bold_italic=texfont(:bolditalic)),
		palette = (color = ColorSchemes.Set1_3, marker = [:circle]),
		Axis = (xlabelsize = 20, ylabelsize = 20, xgridvisible = false, ygridvisible = false, titlesize = 28,
			 xtickalign = 1, ytickalign = 1, xticksmirrored = true, yticksmirrored = true),
		Legend = (labelsize = 12, framevisible = false, ),
		Lines = (linewidth = 3, color = :black, linestyle = :solid, ),
		Scatter = (strokewidth = 1, strokecolor = :black, markersize = 16,)
	)
	set_theme!(My_theme)
end
