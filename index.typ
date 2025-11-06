// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = line(start: (25%,0%), end: (75%,0%))

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): set block(
    fill: luma(230),
    width: 100%,
    inset: 8pt,
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.abs
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == str {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == content {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != str {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    block(below: 0pt, new_title_block) +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black, body_background_color: white) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: body_background_color, width: 100%, inset: 8pt, body))
      }
    )
}



#let article(
  title: none,
  subtitle: none,
  authors: none,
  date: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  lang: "en",
  region: "US",
  font: "libertinus serif",
  fontsize: 11pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: "libertinus serif",
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  sectionnumbering: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)
  if title != none {
    align(center)[#block(inset: 2em)[
      #set par(leading: heading-line-height)
      #if (heading-family != none or heading-weight != "bold" or heading-style != "normal"
           or heading-color != black) {
        set text(font: heading-family, weight: heading-weight, style: heading-style, fill: heading-color)
        text(size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(size: subtitle-size)[#subtitle]
        }
      } else {
        text(weight: "bold", size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(weight: "bold", size: subtitle-size)[#subtitle]
        }
      }
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)
#import "@preview/fontawesome:0.5.0": *

#set page(
  paper: "us-letter",
  margin: (x: 1in,y: 1in,),
  numbering: "1",
)

#show: doc => article(
  title: [CEVE 543 Fall 2025 Lab 7: Bias Correction Implementation],
  subtitle: [Delta method and quantile mapping for temperature bias correction],
  authors: (
    ( name: [Alisa Sannikova],
      affiliation: [],
      email: [] ),
    ),
  date: [2025-10-24],
  fontsize: 11pt,
  sectionnumbering: "1.1.a",
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 1,
  doc,
)

= Background and Goals
<background-and-goals>
Climate models have systematic biases that directly affect impact assessments, arising from coarse spatial resolution, parameterization of sub-grid processes, representation of topography, and errors in simulated circulation patterns. This lab implements two widely-used bias correction approaches: the delta method and quantile-quantile (QQ) mapping. The delta method preserves the climate model's change signal while anchoring absolute values to observations, whereas QQ-mapping corrects the full distribution of values.

Both methods assume stationarity---that the statistical relationship between model and observations remains constant across climate states. This assumption may not hold under significant climate change. We'll explore the strengths and limitations of each method using temperature data for Boston, providing hands-on experience before PS2 Part 1.

= Study Location and Data
<study-location-and-data>
This lab uses temperature data for Boston Logan International Airport (Station ID: USW00014739, 42.3631$""^circle.stroked.tiny$N, 71.0064$""^circle.stroked.tiny$W). Observational data comes from GHCN-Daily (1936-2024) and is provided in the `USW00014739.csv` file in degrees Celsius. Climate model data comes from the GFDL-ESM4 model's 3-hourly near-surface air temperature (`tas`), pre-downloaded from Google Cloud Storage for both historical (1850-2014) and SSP3-7.0 (2015-2100) scenarios. Refer to #link("./labs/Lab-6/index.qmd")[Lab 6] for details on CMIP6 data structure.

== Data Processing Notes
<data-processing-notes>
The GHCN data provides daily average temperature in degrees Celsius. CMIP6 provides 3-hourly instantaneous temperature in Kelvin, which we'll convert to Celsius and aggregate to daily averages. The pre-downloaded NetCDF files (`boston_historical.nc` and `boston_ssp370.nc`) contain 3-hourly surface air temperature for the grid cell nearest to Boston. We aggregate 8 consecutive 3-hour periods to approximate daily averages, ignoring daylight saving time---a simplification typical in many bias correction applications. If you're interested in applying these methods to a different location, the `download_data.jl` script demonstrates how to extract CMIP6 data from Google Cloud Storage.

#block[
#callout(
body: 
[
Before starting the lab, uncomment the `Pkg.instantiate()` line in the first code block and run it to install all required packages. This will take a few minutes the first time. After installation completes, comment the line back out to avoid reinstalling on subsequent runs.

]
, 
title: 
[
Before Starting
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
= Lab Implementation
<lab-implementation>
== Package Setup
<package-setup>
#block[
```julia
using Pkg
lab_dir = dirname(@__FILE__)
Pkg.activate(lab_dir)
# Pkg.instantiate() # uncomment this the first time you run the lab to install packages, then comment it back
```

]
#block[
```julia
using CSV, CairoMakie, DataFrames, Dates, LaTeXStrings, NCDatasets, Statistics, StatsBase, TidierData
ENV["DATAFRAMES_ROWS"] = 6
CairoMakie.activate!()

# Constants for plotting
const MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
const MONTH_LABELS = (1:12, MONTH_NAMES)
```

]
== Task 1: Load and Process Data
<task-1-load-and-process-data>
We begin by loading both observational and climate model data. The observational data from GHCN-Daily spans 1936-2024, while the CMIP6 data requires processing from 3-hourly to daily resolution. We'll use the overlapping historical period (1995-2014) to calibrate bias corrections.

#block[
#callout(
body: 
[
Load the Boston GHCN temperature data from `USW00014739.csv` and filter to years with at least 80% complete data. Then load the pre-downloaded CMIP6 NetCDF files and aggregate the 3-hourly temperature data to daily values for the historical (1995-2014) and near-term future (2020-2040) periods. The helper functions `load_cmip6_data()` and `aggregate_to_daily()` are provided below. Finally, visualize the annual cycle comparing observations vs the historical GCM simulation.

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
#block[
```julia
# Load and clean the observational data
data_path = joinpath(lab_dir, "USW00014739.csv")
df = @chain begin
    CSV.read(data_path, DataFrame)
    @mutate(
        TAVG = ifelse.(ismissing.(TMIN) .| ismissing.(TMAX), missing, (TMIN .+ TMAX) ./ 2),
    )
    @mutate(TAVG = TAVG / 10.0) # Convert to degrees C
    @rename(date = DATE)
    @mutate(year = year(date), month = month(date))
    @select(date, year, month, TAVG)
end

# Filter to years with at least 80% complete data
yearly_counts = @chain df begin
    @group_by(year)
    @summarize(n_obs = sum(!ismissing(TAVG)), n_total = n())
    @mutate(frac_complete = n_obs / n_total)
end

good_years_df = @chain yearly_counts begin
    @filter(frac_complete >= 0.8)
end
good_years = good_years_df.year

df_clean = @chain df begin
    @filter(year in !!good_years)
    dropmissing(:TAVG)
end
```

]
#block[
```julia
# Helper functions for CMIP6 data
"""
Load CMIP6 temperature data from a local NetCDF file and convert times to DateTime.
"""
function load_cmip6_data(file_path::String)
    ds = NCDataset(file_path)
    tas_data = ds["tas"][:]
    time_cf = ds["time"][:]
    close(ds)

    time_data = [DateTime(
        Dates.year(t), Dates.month(t), Dates.day(t),
        Dates.hour(t), Dates.minute(t), Dates.second(t)
    ) for t in time_cf]

    return tas_data, time_data
end

"""
Aggregate 3-hourly temperature data to daily averages.
"""
function aggregate_to_daily(tas_3hr, time_3hr)
    n_3hr_per_day = 8
    n_days = div(length(tas_3hr), n_3hr_per_day)

    daily_temp = Vector{Float64}(undef, n_days)
    daily_dates = Vector{Date}(undef, n_days)

    for i in 1:n_days
        idx_start = (i - 1) * n_3hr_per_day + 1
        idx_end = i * n_3hr_per_day
        daily_vals = collect(skipmissing(tas_3hr[idx_start:idx_end]))
        daily_temp[i] = isempty(daily_vals) ? NaN : mean(daily_vals)
        daily_dates[i] = Date(time_3hr[idx_start])
    end

    daily_temp_c = daily_temp .- 273.15  # Convert K to C
    return daily_temp_c, daily_dates
end

# Load and process CMIP6 data
hist_file = joinpath(lab_dir, "boston_historical.nc")
ssp370_file = joinpath(lab_dir, "boston_ssp370.nc")

tas_hist_3hr, time_hist_3hr = load_cmip6_data(hist_file)
tas_ssp370_3hr, time_ssp370_3hr = load_cmip6_data(ssp370_file)

# Process historical period (1995-2014)
hist_start = DateTime(1995, 1, 1)
hist_end = DateTime(2014, 12, 31, 23, 59, 59)
hist_idx = (time_hist_3hr .>= hist_start) .& (time_hist_3hr .<= hist_end)
tas_hist_daily, dates_hist_daily = aggregate_to_daily(tas_hist_3hr[hist_idx], time_hist_3hr[hist_idx])

# Process near-term period (2020-2040)
near_start = DateTime(2020, 1, 1)
near_end = DateTime(2040, 12, 31, 23, 59, 59)
near_idx = (time_ssp370_3hr .>= near_start) .& (time_ssp370_3hr .<= near_end)
tas_ssp370_near_daily, dates_ssp370_near_daily = aggregate_to_daily(tas_ssp370_3hr[near_idx], time_ssp370_3hr[near_idx])

# Create DataFrames
df_gcm_hist = DataFrame(
    date=dates_hist_daily, temp=tas_hist_daily,
    year=year.(dates_hist_daily), month=month.(dates_hist_daily)
)

df_ssp370_near = DataFrame(
    date=dates_ssp370_near_daily, temp=tas_ssp370_near_daily,
    year=year.(dates_ssp370_near_daily), month=month.(dates_ssp370_near_daily)
)

df_obs_hist = @chain df_clean begin
    @filter(year >= 1995 && year <= 2014)
end

# Compute monthly climatologies
obs_monthly = @chain df_obs_hist begin
    @group_by(month)
    @summarize(mean_temp = mean(TAVG))
end

gcm_monthly = @chain df_gcm_hist begin
    @group_by(month)
    @summarize(mean_temp = mean(temp))
end
```

]
```julia
let
    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel="Month",
        ylabel=L"Temperature ($^\circ$C)",
        title="Annual Cycle: Observations vs GCM Historical (1995-2014)",
        xticks=MONTH_LABELS)
    lines!(ax, obs_monthly.month, obs_monthly.mean_temp, linewidth=2, color=:steelblue, label="Observations")
    lines!(ax, gcm_monthly.month, gcm_monthly.mean_temp, linewidth=2, color=:coral, label="GCM Historical")
    axislegend(ax, position=:lt)
    fig
end
```

#figure([
#box(image("index_files/figure-typst/fig-gcm-obs-comparison-output-1.png", height: 5in, width: 7in))
], caption: figure.caption(
position: bottom, 
[
Annual cycle comparison between GHCN observations and GFDL-ESM4 historical simulation for Boston (1995-2014). The GCM shows a warm bias across most months.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-gcm-obs-comparison>


== Task 2: Implement Delta Method
<task-2-implement-delta-method>
The delta method corrects the mean bias while preserving the climate model's projected change signal. We calculate a monthly bias correction based on the historical period (1995-2014) and apply it to future projections.

#block[
#callout(
body: 
[
Implement the additive delta method for temperature bias correction. For each calendar month $m$, calculate the mean bias: $Delta_m = macron(T)_(upright("hist") \, m)^(upright("GCM")) - macron(T)_(upright("hist") \, m)^(upright("obs"))$. Then apply the correction to future values: $T_(upright("fut"))^(upright("corr")) (d \, m \, y) = T_(upright("fut"))^(upright("GCM")) (d \, m \, y) - Delta_m$.

Follow these steps:

+ Calculate the monthly mean bias by grouping both `df_gcm_hist` and `df_obs_hist` by month, computing their means, joining them, and computing the difference.
+ Create a bar plot visualizing the monthly bias.
+ Write a function `apply_delta_method(gcm_temps, gcm_dates, monthly_bias_df)` that applies the bias correction to a vector of temperatures and dates.
+ Apply your function to the near-term data and add a new column `temp_delta` to `df_ssp370_near`.
+ Create monthly climatologies for both raw and delta-corrected near-term data.
+ Visualize the annual cycle showing historical observations, raw GCM near-term, and delta-corrected near-term temperatures.

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
```julia
#create a dataframe
monthly_bias = @chain leftjoin(gcm_monthly, obs_monthly; on = :month, makeunique = true) begin
    rename!(:mean_temp => :gcm_mean, :mean_temp_1 => :obs_mean)
    @mutate(bias = gcm_mean .- obs_mean)
end
```

#table(
  columns: 5,
  align: (right,right,right,right,right,),
  table.header(table.cell(align: right)[#set text(weight: "bold"); Row], table.cell(align: left)[month], table.cell(align: left)[gcm\_mean], table.cell(align: left)[obs\_mean], table.cell(align: left)[bias],
    table.cell(align: right)[#set text(weight: "bold"); ], table.cell(align: left)[Int64], table.cell(align: left)[Float64], table.cell(align: left)[Float64?], table.cell(align: left)[Float64],),
  table.hline(),
  table.cell(align: right)[#set text(weight: "bold"); 1], table.cell(align: right)[1], table.cell(align: right)[-3.18274], table.cell(align: right)[-1.09024], table.cell(align: right)[-2.0925],
  table.cell(align: right)[#set text(weight: "bold"); 2], table.cell(align: right)[2], table.cell(align: right)[-1.70476], table.cell(align: right)[0.0218584], table.cell(align: right)[-1.72662],
  table.cell(align: right)[#set text(weight: "bold"); 3], table.cell(align: right)[3], table.cell(align: right)[2.00539], table.cell(align: right)[3.78944], table.cell(align: right)[-1.78405],
  table.cell(align: right)[⋮], table.cell(align: right)[⋮], table.cell(align: right)[⋮], table.cell(align: right)[⋮], table.cell(align: right)[⋮],
  table.cell(align: right)[#set text(weight: "bold"); 10], table.cell(align: right)[10], table.cell(align: right)[11.555], table.cell(align: right)[12.6931], table.cell(align: right)[-1.1381],
  table.cell(align: right)[#set text(weight: "bold"); 11], table.cell(align: right)[11], table.cell(align: right)[5.14819], table.cell(align: right)[7.06367], table.cell(align: right)[-1.91548],
  table.cell(align: right)[#set text(weight: "bold"); 12], table.cell(align: right)[12], table.cell(align: right)[0.52005], table.cell(align: right)[2.04863], table.cell(align: right)[-1.52858],
)
```julia
# Create a new figure
fig = Figure(resolution = (800, 400))

# Add an axis
ax = Axis(fig[1, 1];
    xticks = (1:12, ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]),
    xlabel = "Month",
    ylabel = "Temperature Bias (°C)",
    title  = "Monthly Mean Temperature Bias"
)

# Plot bars
barplot!(ax, monthly_bias.month, monthly_bias.bias, color = :steelblue)

# Add a horizontal line at 0
hlines!(ax, 0, color = :black, linewidth = 1.5)

# Show the figure (last line ensures it renders)
fig
```

#block[
```
┌ Warning: Found `resolution` in the theme when creating a `Scene`. The `resolution` keyword for `Scene`s and `Figure`s has been deprecated. Use `Figure(; size = ...` or `Scene(; size = ...)` instead, which better reflects that this is a unitless size and not a pixel resolution. The key could also come from `set_theme!` calls or related theming functions.
└ @ Makie ~/.julia/packages/Makie/4JW9B/src/scenes.jl:264
```

]
#box(image("index_files/figure-typst/cell-8-output-2.png", height: 4.16667in, width: 8.33333in))

```julia
function apply_delta_method(gcm_temps, gcm_dates, monthly_bias_df)
    # Create vector to store corrected temperatures
    corrected_temps = similar(gcm_temps)

    # Loop through each temperature and date
    for i in eachindex(gcm_temps)
        m = month(gcm_dates[i])  # extract month number (1–12)
        # find the bias for that month
        bias_val = monthly_bias_df.bias[monthly_bias_df.month .== m][1]
        # subtract bias from model temperature
        corrected_temps[i] = gcm_temps[i] - bias_val
    end

    return corrected_temps
end
```

```
apply_delta_method (generic function with 1 method)
```

```julia
df_ssp370_near.temp_delta = apply_delta_method(
    df_ssp370_near.temp,
    df_ssp370_near.date,
    monthly_bias
)
```

```
7665-element Vector{Float64}:
 -5.293100300450412
 -5.623239460606662
 -0.9786776930285366
  2.5020046800183384
  5.411062297205838
  0.4519863694714634
 -3.4730624586535366
 -6.282297077794162
  1.1962185472058384
 -2.8390292555285366
  ⋮
  2.642983674080141
  0.09229397681451634
 -0.6736667165448587
  0.44959378150201634
 -1.1007602224042337
 -1.8960483083417337
  1.2649624338457663
 -1.7724215993573587
 -3.997549773185484
```

```julia
near_monthly_raw = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp_raw = mean(temp))
end

near_monthly_delta = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp_delta = mean(temp_delta))
end
```

#table(
  columns: 3,
  align: (right,right,right,),
  table.header(table.cell(align: right)[#set text(weight: "bold"); Row], table.cell(align: left)[month], table.cell(align: left)[mean\_temp\_delta],
    table.cell(align: right)[#set text(weight: "bold"); ], table.cell(align: left)[Int64], table.cell(align: left)[Float64],),
  table.hline(),
  table.cell(align: right)[#set text(weight: "bold"); 1], table.cell(align: right)[1], table.cell(align: right)[0.103505],
  table.cell(align: right)[#set text(weight: "bold"); 2], table.cell(align: right)[2], table.cell(align: right)[1.70296],
  table.cell(align: right)[#set text(weight: "bold"); 3], table.cell(align: right)[3], table.cell(align: right)[4.62097],
  table.cell(align: right)[⋮], table.cell(align: right)[⋮], table.cell(align: right)[⋮],
  table.cell(align: right)[#set text(weight: "bold"); 10], table.cell(align: right)[10], table.cell(align: right)[13.4868],
  table.cell(align: right)[#set text(weight: "bold"); 11], table.cell(align: right)[11], table.cell(align: right)[7.91067],
  table.cell(align: right)[#set text(weight: "bold"); 12], table.cell(align: right)[12], table.cell(align: right)[2.51416],
)
```julia
fig = Figure(resolution = (900, 420))
ax = Axis(fig[1, 1];
    xticks = MONTH_LABELS,  # uses the provided constant
    xlabel = "Month",
    ylabel = "Temperature (°C)",
    title = "Delta Method: Annual Cycle for SSP3-7.0 Near-Term (2020–2040)"
)

# Historical Obs (blue)
lines!(ax, obs_monthly.month, obs_monthly.mean_temp,
    color = :blue, label = "Historical Obs", linewidth = 2)

# GCM Raw (coral/red)
lines!(ax, near_monthly_raw.month, near_monthly_raw.mean_temp_raw,
    color = :coral, label = "GCM Raw", linewidth = 2)

# Delta Corrected (green dashed)
lines!(ax, near_monthly_delta.month, near_monthly_delta.mean_temp_delta,
    color = :green, linestyle = :dash, label = "Delta Corrected", linewidth = 2)

axislegend(ax, position = :rb)  # bottom-right legend
fig
```

#block[
```
┌ Warning: Found `resolution` in the theme when creating a `Scene`. The `resolution` keyword for `Scene`s and `Figure`s has been deprecated. Use `Figure(; size = ...` or `Scene(; size = ...)` instead, which better reflects that this is a unitless size and not a pixel resolution. The key could also come from `set_theme!` calls or related theming functions.
└ @ Makie ~/.julia/packages/Makie/4JW9B/src/scenes.jl:264
```

]
#box(image("index_files/figure-typst/cell-12-output-2.png", height: 4.375in, width: 9.375in))

== Task 3: Implement Quantile-Quantile Mapping
<task-3-implement-quantile-quantile-mapping>
Unlike the delta method, QQ-mapping transforms the entire probability distribution of model output to match observations. For each value in the future model output, we find its percentile in the historical model distribution, then map it to the same percentile in the historical observed distribution.

#block[
#callout(
body: 
[
Implement QQ-mapping to correct the full distribution of temperature values.

Follow these steps:

+ Extract the observed and GCM historical temperatures as vectors (`obs_hist_temps` and `gcm_hist_temps`).
+ Use `ecdf()` from StatsBase.jl to fit empirical cumulative distribution functions to both datasets.
+ Write a function `apply_qqmap_empirical(gcm_temps, ecdf_gcm, obs_hist_temps)` that:
  - For each temperature value, finds its percentile in the GCM CDF
  - Clamps the percentile to the range \[0.001, 0.999\] to avoid extrapolation
  - Maps to the same percentile in the observed distribution using `quantile()`
+ Apply your function to the near-term data and add a new column `temp_qqmap` to `df_ssp370_near`.
+ Create a QQ-plot comparing observed and GCM quantiles for the historical period, showing the 1:1 line.

#block[
#callout(
body: 
[
- Use `ecdf_gcm(value)` to get the percentile of a value in the GCM distribution
- Use `quantile(obs_hist_temps, p)` to get the value at percentile `p` in the observed distribution
- The `clamp(x, low, high)` function constrains `x` to the range \[low, high\]

]
, 
title: 
[
Hints
]
, 
background_color: 
rgb("#ccf1e3")
, 
icon_color: 
rgb("#00A047")
, 
icon: 
fa-lightbulb()
, 
body_background_color: 
white
)
]
]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
```julia
# Extract temperature vectors
obs_hist_temps = Float64.(df_obs_hist.TAVG)
gcm_hist_temps = Float64.(df_gcm_hist.temp)

ecdf_obs = ecdf(obs_hist_temps)
ecdf_gcm = ecdf(gcm_hist_temps)
```

```
ECDF{Vector{Float64}, Weights{Float64, Float64, Vector{Float64}}}([-17.170874023437477, -16.842443847656227, -16.635809326171852, -16.161810302734352, -16.041326904296852, -15.905065917968727, -15.873297119140602, -15.437841796874977, -15.356787109374977, -15.301611328124977  …  26.878961181640648, 26.895562744140648, 26.924462890625023, 26.955651855468773, 27.055413818359398, 27.076562500000023, 27.093164062500023, 27.192712402343773, 27.408654785156273, 27.778253173828148], Float64[])
```

```julia
function apply_qqmap_empirical(gcm_temps, ecdf_gcm, obs_hist_temps)
    # result vector (Float64)
    corrected_temps = similar(gcm_temps, Float64)

    for i in eachindex(gcm_temps)
        # percentile of the model value in the model's historical distribution
        p = ecdf_gcm(gcm_temps[i])

        # clamp to avoid exactly 0 or 1 (safer interpolation at the tails)
        p = clamp(p, 0.001, 0.999)

        # map to the same percentile in the observed historical distribution
        corrected_temps[i] = quantile(obs_hist_temps, p)
    end

    return corrected_temps
end
```

```
apply_qqmap_empirical (generic function with 1 method)
```

```julia
# Apply QQ-mapping correction to the future GCM data
df_ssp370_near.temp_qqmap = apply_qqmap_empirical(
    df_ssp370_near.temp,
    ecdf_gcm,
    obs_hist_temps
)
```

```
7665-element Vector{Float64}:
 -5.75
 -6.1
 -1.95
  2.5
  5.8
 -0.3
 -4.15
 -6.65
  0.8
 -3.65
  ⋮
  3.35
  0.0
 -1.1
  0.55
 -1.65
 -2.25
  1.65
 -2.2
 -4.15
```

```julia
# Define percentiles
percentiles = 0.01:0.01:0.99

# Compute quantiles for observed and model historical data
obs_quantiles = quantile(obs_hist_temps, percentiles)
gcm_hist_quantiles = quantile(gcm_hist_temps, percentiles)

# Create the figure
fig = Figure(resolution = (600, 600))
ax = Axis(fig[1, 1],
    xlabel = "Observed Quantiles (°C)",
    ylabel = "GCM Quantiles (°C)",
    title = "QQ Plot: Historical GCM vs Observed",
    aspect = DataAspect()   # makes it square
)

# Scatter plot of quantiles
scatter!(ax, obs_quantiles, gcm_hist_quantiles, color = :steelblue, markersize = 6)

# 1:1 reference line (perfect match)
lines!(ax, [-20, 40], [-20, 40], color = :black, linestyle = :dash)

fig
```

#block[
```
┌ Warning: Found `resolution` in the theme when creating a `Scene`. The `resolution` keyword for `Scene`s and `Figure`s has been deprecated. Use `Figure(; size = ...` or `Scene(; size = ...)` instead, which better reflects that this is a unitless size and not a pixel resolution. The key could also come from `set_theme!` calls or related theming functions.
└ @ Makie ~/.julia/packages/Makie/4JW9B/src/scenes.jl:264
```

]
#box(image("index_files/figure-typst/cell-16-output-2.png", height: 6.25in, width: 6.25in))

== Task 4: Compare Methods
<task-4-compare-methods>
Now we compare the delta method and QQ-mapping approaches to understand their strengths and limitations.

#block[
#callout(
body: 
[
Create a single figure showing monthly mean temperature for the near-term period (2020-2040) using all four approaches:

+ Compute the monthly mean for QQ-mapped temperatures (`near_monthly_qqmap`)
+ Create a figure with 4 lines:
  - Historical observations (`obs_monthly`)
  - Raw GCM near-term (`near_monthly_raw`)
  - Delta-corrected near-term (`near_monthly_delta`)
  - QQ-mapped near-term (`near_monthly_qqmap`)

After creating the plot, examine it carefully and consider:

- How do the two correction methods differ?
- What does the QQ-mapped line reveal about this method's treatment of the warming signal?

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
```julia
near_monthly_qqmap = @chain df_ssp370_near begin
    @group_by(month)
    @summarize(mean_temp_qqmap = mean(temp_qqmap))
end
```

#table(
  columns: 3,
  align: (right,right,right,),
  table.header(table.cell(align: right)[#set text(weight: "bold"); Row], table.cell(align: left)[month], table.cell(align: left)[mean\_temp\_qqmap],
    table.cell(align: right)[#set text(weight: "bold"); ], table.cell(align: left)[Int64], table.cell(align: left)[Float64],),
  table.hline(),
  table.cell(align: right)[#set text(weight: "bold"); 1], table.cell(align: right)[1], table.cell(align: right)[-0.374928],
  table.cell(align: right)[#set text(weight: "bold"); 2], table.cell(align: right)[2], table.cell(align: right)[1.75688],
  table.cell(align: right)[#set text(weight: "bold"); 3], table.cell(align: right)[3], table.cell(align: right)[4.83506],
  table.cell(align: right)[⋮], table.cell(align: right)[⋮], table.cell(align: right)[⋮],
  table.cell(align: right)[#set text(weight: "bold"); 10], table.cell(align: right)[10], table.cell(align: right)[14.1417],
  table.cell(align: right)[#set text(weight: "bold"); 11], table.cell(align: right)[11], table.cell(align: right)[8.04862],
  table.cell(align: right)[#set text(weight: "bold"); 12], table.cell(align: right)[12], table.cell(align: right)[2.77888],
)
```julia
# Create figure and axis
fig = Figure(resolution = (900, 420))
ax = Axis(fig[1, 1],
    xlabel = "Month",
    ylabel = "Temperature (°C)",
    title = "Annual Cycle Comparison: All Methods (Near-term 2020–2040)",
    xticks = MONTH_LABELS,   # use provided month names
)

# Plot lines for each dataset
lines!(ax, obs_monthly.month, obs_monthly.mean_temp,
    color = :black, linewidth = 2.5, label = "Historical Obs")

lines!(ax, near_monthly_raw.month, near_monthly_raw.mean_temp_raw,
    color = :coral, linewidth = 2, label = "Raw GCM")

lines!(ax, near_monthly_delta.month, near_monthly_delta.mean_temp_delta,
    color = :green, linewidth = 2, linestyle = :dash, label = "Delta")

lines!(ax, near_monthly_qqmap.month, near_monthly_qqmap.mean_temp_qqmap,
    color = :purple, linewidth = 2, linestyle = :dot, label = "QQ-map")

# Add legend
axislegend(ax, position = :lt)

fig
```

#block[
```
┌ Warning: Found `resolution` in the theme when creating a `Scene`. The `resolution` keyword for `Scene`s and `Figure`s has been deprecated. Use `Figure(; size = ...` or `Scene(; size = ...)` instead, which better reflects that this is a unitless size and not a pixel resolution. The key could also come from `set_theme!` calls or related theming functions.
└ @ Makie ~/.julia/packages/Makie/4JW9B/src/scenes.jl:264
```

]
#box(image("index_files/figure-typst/cell-18-output-2.png", height: 4.375in, width: 9.375in))

== Task 5: Reflection
<task-5-reflection>
Finally, we reflect on the assumptions, limitations, and appropriate use cases for each bias correction method.

#block[
#callout(
body: 
[
Write brief responses (2-3 sentences each) to the following questions:

+ #strong[Method selection:] If you were providing climate data to support a decision about urban heat management in Boston, which bias correction method would you recommend and why?

+ #strong[Appropriateness of QQ-mapping:] In the #cite(<ines_biascorrection:2006>, form: "prose") paper we discussed, QQ-mapping was used to correct both the #emph[frequency] and #emph[intensity] of rainfall for crop modeling. For temperature data, does it make sense to correct the full distribution? Why or why not? When might the delta method be more appropriate than QQ-mapping?

]
, 
title: 
[
Instructions
]
, 
background_color: 
rgb("#f7dddc")
, 
icon_color: 
rgb("#CC1914")
, 
icon: 
fa-exclamation()
, 
body_background_color: 
white
)
]
=== Method Selection for Urban Heat Management
<method-selection-for-urban-heat-management>
The models showed fairly similar results. However, I believe that the delta method should be used, as it better accounts for long-term climate change shifts. In contrast, QQ-mapping could incorrectly attribute a fixed percentile relationship to a distribution that is constantly changing. Therefore, for long records such as 40 years, the delta method is more appropriate.

=== Appropriateness of QQ-Mapping for Temperature
<appropriateness-of-qq-mapping-for-temperature>
The QQ method better represents extreme events, allowing short-term (e.g., weekly) variations to be modeled more accurately.

= References
<references>




