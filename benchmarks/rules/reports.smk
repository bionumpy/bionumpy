import plotly.express as px

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


def get_runtime_from_benchmark_report(file_name):
    tool_name = file_name.split("/")[2]
    time = float(open(file_name).readlines()[1].split()[0])
    return tool_name, time


def get_report_name(wildcards):
    return config["analysis"][wildcards.analysis]["name"]


rule make_runtime_report:
    input:
        get_analysis_benchmark_reports
    output:
        data="results/reports/{analysis}/{dataset_size}.data",
        figure="results/reports/{analysis}/{dataset_size}.html",
        png="results/reports/{analysis}/{dataset_size}.png"
    params:
        report_name=get_report_name
    run:
        print(input)
        from collections import defaultdict
        runtimes = defaultdict(float)
        for benchmark_file in input:
            method_name, time = get_runtime_from_benchmark_report(benchmark_file)
            pretty_name = config["method_names"][method_name]
            runtimes[pretty_name] += time

        print(runtimes)
        with open(output.data, "w") as f:
            f.write("\n".join(["\t".join(map(str, result)) for result in runtimes.items()]))

        x_axis = runtimes.keys()  # [r for result in results]
        y_axis = runtimes.values()  # [result[1] for result in results]
        fig = px.bar(x=x_axis, y=y_axis, title=params.report_name,
            labels={"x": "Method",
                    "y": "Time in seconds"
                    }
        )
        fig.write_html(output.figure)
        fig.write_image(output.png, scale=4)


rule main_report:
    input:
        ["results/reports/" + analysis + "/{run_size}.png" for analysis in config["analysis"]]
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
            if i % 2 == 0:
                out += "<br>";
            out += "<img src='" + image + "' style='width: 50%; height: auto'>"

        with open(output.html, "w") as f:
            f.write(out)
