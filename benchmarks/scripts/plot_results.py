import plotly.express as px
from collections import defaultdict

def get_runtime_from_benchmark_report(file_name):
    tool_name = file_name.split("/")[2]
    time = float(open(file_name).readlines()[1].split()[0])
    return tool_name, time


runtimes = defaultdict(float)
for benchmark_file in snakemake.input:
    method_name, time = get_runtime_from_benchmark_report(benchmark_file)
    pretty_name = snakemake.config["method_names"][method_name]
    runtimes[pretty_name] += time

with open(snakemake.output.data, "w") as f:
    f.write("\n".join(["\t".join(map(str, result)) for result in runtimes.items()]))

x_axis = runtimes.keys()  # [r for result in results]
y_axis = runtimes.values()  # [result[1] for result in results]
fig = px.bar(x=x_axis, y=y_axis, title=snakemake.params.report_name,
    labels={"x": "Method",
            "y": "Time in seconds"
            },
    template="seaborn",
)
fig.update_layout(
    font={
        "size": 19
    }
)
fig.write_html(snakemake.output.figure)
fig.write_json(snakemake.output.json)
fig.write_image(snakemake.output.png, scale=2)


