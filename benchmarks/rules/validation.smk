import itertools


def get_result_files_for_two_methods(wildcards):
    analysis_config = config["analysis"][wildcards.analysis]
    file_ending = analysis_config["output_file_ending"]

    return ["results/" + method + "/" + wildcards.analysis + "/" + wildcards.run_details + file_ending
            for method in (wildcards.method1, wildcards.method2)]


# checking that two methods give the same results
rule assert_method_output_is_equal:
    input:
        get_result_files_for_two_methods
    output:
        "results/validation/{analysis}/{method1}-vs-{method2}.{run_details}.txt"
    conda:
        "../envs/diff.yml"
    shell:
        # diff return nonzero exit code on diff
        """
        diff --ignore-space-change {input}
        touch {output}
        """


def get_validation_files_for_analysis(wildcards):
    tools = config["analysis"][wildcards.analysis]["tools"]
    tool_pairs = itertools.combinations(tools, 2)
    out = ["results/validation/" + wildcards.analysis + "/" + tool1 + "-vs-" + tool2 + "." + wildcards.run_details + ".txt"
           for tool1, tool2 in tool_pairs]
    print(out)
    return out


# validates that all methods give same result for an analysis
rule single_validation:
    input:
        get_validation_files_for_analysis
    output:
        "results/validation/{analysis}/{run_details}.validation"
    shell:
        """
        cat {input} > {output}
        """


def get_validation_reports(wildcards):
    # finds all types of analysis that should be validated
    out = []
    for analysis in config["analysis"]:
        if config["analysis"][analysis]["validate_equal"]:
            run_details = config["analysis"][analysis]["runs"][wildcards.dataset_size]
            out.append("results/validation/" + analysis + "/" + run_details + ".validation")

    print(out)
    return out


rule validation_report_all_analysis:
    input:
        get_validation_reports
    output:
        "validation_report_{dataset_size}.md"
    shell:
        """
        cat {input} > {output}
        echo 'passed\n' >> {output}
        """
