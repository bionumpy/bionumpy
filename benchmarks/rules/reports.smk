
def get_analysis_benchmark_reports(wildcards):
    analysis_type = wildcards.analysis
    analysis_config = config["analysis"][analysis_type]
    assert wildcards.dataset_size in ["small", "big"], wildcards.dataset_size

    dataset_size = analysis_config["runs"][wildcards.dataset_size]
    tools = analysis_config["tools"]
    return [f"benchmarks/" + wildcards.analysis + "/" + tool + "/" + dataset_size + ".txt" for tool in tools]


def get_runtime_from_benchmark_report(file_name):
    tool_name = file_name.split("/")[2]
    time = float(open(file_name).readlines()[1].split()[0])
    return tool_name, time


rule make_runtime_report:
    input:
        get_analysis_benchmark_reports
    output:
        data="results/reports/{analysis}/{dataset_size}.data",
        figure="results/reports/{analysis}/{dataset_size}.html",
        png="results/reports/{analysis}/{dataset_size}.png"
    run:
        results = [get_runtime_from_benchmark_report(i) for i in input]
        with open(output.data, "w") as f:
            f.write("\n".join(["\t".join(map(str, result)) for result in results]))

        x_axis = [result[0] for result in results]
        y_axis = [result[1] for result in results]
        fig = px.bar(x=x_axis, y=y_axis, title=wildcards.analysis)
        fig.write_html(output.figure)
        fig.write_image(output.png, scale=4)

