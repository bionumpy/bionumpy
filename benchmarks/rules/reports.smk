
def get_analysis_benchmark_reports(wildcards):
    analysis_type = wildcards.analysis
    analysis_config = config["analysis"][analysis_type]
    assert wildcards.dataset_size in ["small", "big"], wildcards.dataset_size

    dataset_size = analysis_config["runs"][wildcards.dataset_size]
    tools = analysis_config["tools"]
    
    files = []
    for tool in tools:
        benchmark_prefixes = [""]
        if "benchmark_report_prefixes" in config["analysis"][analysis_type]:
            if tool in config["analysis"][analysis_type]["benchmark_report_prefixes"]:
                benchmark_prefixes = config["analysis"][analysis_type]["benchmark_report_prefixes"][tool]

        for prefix in benchmark_prefixes:
            files.append("benchmarks/" + wildcards.analysis + "/" + tool + "/" + prefix + dataset_size + ".txt")

    return files



def get_report_name(wildcards):
    return config["analysis"][wildcards.analysis]["name"]


rule make_runtime_report:
    input:
        get_analysis_benchmark_reports
    output:
        data="results/reports/{analysis}/{dataset_size}.data",
        figure="results/reports/{analysis}/{dataset_size}.html",
        png="results/reports/{analysis}/{dataset_size}.png",
        json="results/reports/{analysis}/{dataset_size}.json"
    params:
        report_name=get_report_name
    script:
        "../scripts/plot_results.py"


rule main_report:
    input:
        ["results/reports/" + analysis + "/{run_size}.png" for analysis in config["analysis"]],
        'validation_report_{run_size}.md'
    output:
        md="report_{run_size}.md",
        html="report_{run_size}.html"
    run:
        # markdown
        out = "# Benchmark report \n\n"
        out += "\n\n".join("![](" + image + ")" for image in input)
        with open(output.md, "w") as f:
            f.write(out)

        out = ""
        for i, image in enumerate(input):
            if i % 3 == 0:
                out += "<br>";
            out += "<img src='" + image + "' style='width: 33%; height: auto'>"

        with open(output.html, "w") as f:
            f.write(out)

rule main_report_as_png:
    input:
        ["results/reports/" + analysis + "/{run_size}.json" for analysis in config["analysis"]]
    output:
        png="report_{run_size}.png",
    run:
        from plotly.subplots import  make_subplots
        import plotly.io as pio

        n_rows = len(input) // 3

        subfigures = []
        titles = []
        for i, image in enumerate(input):
            subfig = pio.from_json(open(image).read())
            analysis = image.split("/")[2]
            print(analysis)
            title = config["analysis"][analysis]["name"]
            titles.append(title)
            print(type(subfig))
            subfigures.append(subfig.data[0])

        fig = make_subplots(rows=n_rows, cols=3, subplot_titles=titles, y_title="Time in seconds")

        for i, subfig in enumerate(subfigures):
            col = i % 3 + 1
            row = i // 3 + 1
            print(col, row)
            fig.add_trace(subfig, row=row, col=col)

        fig.update_layout(
            font={
                "size": 12
            }
        )

        fig.update_layout(height=1200,width=1500)
        fig.write_image(output[0])
        fig.show()
