def draw_log_dist_plot(df, cell_no, cell_no_ext, cell_no_mrg, options):
	import pandas as pd
	import matplotlib.pyplot as plt
	from matplotlib import colors as mcolors
	import seaborn as sns
	import numpy as np

	ax_x     = df.loc[df['idx'] == cell_no,     "r_log10_idx"].values[0]
	ax_x_ext = df.loc[df['idx'] == cell_no_ext, "r_log10_idx"].values[0]
	ax_x_mrg = df.loc[df['idx'] == cell_no_mrg, "r_log10_idx"].values[0]

	plot_df = pd.concat([pd.DataFrame({'r_log10_idx': df['r_log10_idx'],     \
                                           'value':       df['med_log10_slope'], \
                                           'class':       "med_log10_slope"}),   \
                             pd.DataFrame({'r_log10_idx': df['r_log10_idx'],     \
                                           'value':       df['log10_UMI'],       \
                                           'class':       "log10_UMI"})]).dropna().reset_index().drop_duplicates()

	print("Drawing figure ...", end = "", flush = True)

	sns.set(rc = {"figure.figsize": (12, 7)})
	sns.set_theme(style = "white")

	palette = sns.color_palette(["#b0b0b0", "#000000"])
	enmax_palette = ["#3cc800", "#0100ff"]
	color_codes_wanted = ['my_green3', 'my_blue']
	cdict = dict(zip(color_codes_wanted, [mcolors.to_rgba(c) for c in enmax_palette]))
	mcolors.get_named_colors_mapping().update(cdict)

	sns.lineplot(data = plot_df, x = "r_log10_idx", y = "value", hue = "class", palette = palette).set(xlabel = 'log10(idx)', ylabel = 'log10(UMI number)')

	plt.axvline(x = ax_x,     linestyle = '--', color = 'my_green3')
	plt.axvline(x = ax_x_ext, linestyle = '--', color = 'my_blue')
	plt.axvline(x = ax_x_mrg, linestyle = '-',  color = 'red')

	plt.legend(loc = 'upper right')
	plt.savefig(options.CB_log10_dist_o, format = 'png', dpi = 300)

	print("Done", flush = True)

	return

