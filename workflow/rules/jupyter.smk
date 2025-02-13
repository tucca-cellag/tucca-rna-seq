rule hello:
    output:
        "test.txt",
    conda:
        "../envs/hello.yaml"
    notebook:
        "../notebooks/hello.py.ipynb"
