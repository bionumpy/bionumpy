import numpy as np
import pandas as pd
import plotly.express as px
from collections import defaultdict


def get_runtime_from_benchmark_report(file_name):
    tool_name = file_name.split("/")[2]
    time = float(open(file_name).readlines()[1].split()[0])
    all_times = []
    with open(file_name) as f:
        f.readline()  # skip header
        for line in f:
            time = float(line.split()[0])
            all_times.append(time)
    all_times = np.array(all_times)
    median_time = np.median(all_times)
    print(f"Found median time {median_time} for {tool_name} from {all_times}")
    return tool_name, median_time, all_times


runtimes = defaultdict(float)
all_runtimes = {}
# there might be multiple benchmark times for a single method (if multiple rules are used)
# simply add together the times from each rule
for benchmark_file in snakemake.input:
    method_name, time, all_times = get_runtime_from_benchmark_report(benchmark_file)
    pretty_name = snakemake.config["method_names"][method_name]
    runtimes[pretty_name] += time
    if pretty_name in all_runtimes:
        all_runtimes[pretty_name] += all_times
    else:
        all_runtimes[pretty_name] = all_times


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
# boxplot with all_times
# make a dataframe with method and times
method_names = []
times = []
for method_name, all_times in all_runtimes.items():
    method_names.extend([method_name] * len(all_times))
    times.extend(all_times)

df = pd.DataFrame({"Method": method_names, "Time": times})
fig2 = px.box(df, x="Method", y="Time", title=snakemake.params.report_name,
    labels={"x": "Method",
            "y": "Time in seconds"
            },
    template="seaborn",
)

fig.write_html(snakemake.output.figure)
fig.write_json(snakemake.output.json)
fig.write_image(snakemake.output.png, scale=2)
fig2.write_image(snakemake.output.boxplot, scale=2)
fig2.write_json(snakemake.output.boxplot_json)


