import csv
import matplotlib.pyplot as plt

from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait


chrome_option = webdriver.ChromeOptions()
chrome_option.add_argument(argument='--headless')
browser = webdriver.Chrome(options=chrome_option)   # options=chrome_option
wait = WebDriverWait(browser, 3)  # 最长等待时间为10S
# wait_v = WebDriverWait(browser, 20)
browser.maximize_window()


def load_gene_type():
    gene_name_type_file = open("./mapping_data/combined.gene.bed+2")
    name_type = {}
    for line in gene_name_type_file.readlines():
        info = line.split()[6:]
        g_name = info[0]
        g_type = info[1]
        if g_name not in name_type:
            name_type[g_name] = g_type
    return name_type


def load_mapped_genes():
    gene_file = open("./gene_counts_with_type.csv")
    genes = []
    for line in gene_file.readlines():
        info = line.strip().split(",")
        genes.append(info)
    return genes


def crawler(g_name):
    base_url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    url = base_url + g_name
    browser.get(url)
    print("Searching:", g_name, end="...")
    try:
        g_type = wait.until(EC.presence_of_element_located(
            (By.CSS_SELECTOR, "#aliases_descriptions > div.row "
                              "> div.col-xs-12.col-md-9 > div:nth-child(2) > div > div"))).text
    except TimeoutException:
        try:
            temp = wait.until(EC.presence_of_element_located(
                (By.CSS_SELECTOR, "#cardPage > div:nth-child(1) > div > div > div.gc-section-header.gc-first "
                                  "> div.col-xs-8.col-md-10.col-sm-9.gene-symbol-description.row-no-padding "
                                  "> div > span.gc-category"))).text
            if temp == "Protein Coding":
                g_type = "protein_coding"
            elif temp == "Pseudogene":
                g_type = "pseudogene"
            else:
                print("----------------------------------", temp)
                g_type = temp

        except TimeoutException:
            g_type = "unknow_type"
    print("\t\ttype: ", g_type)
    return g_type


def mapping():
    """
        Available mode: sum, 1, 2, 3
    """
    name_type = load_gene_type()

    tsv_file = open("./mapping_data/gene_counts.tsv")
    my_genes = []
    gene_type_count = {}

    with open("gene_counts_with_type.csv", "w", newline="") as result_file:
        csv_writer = csv.writer(result_file, dialect="excel")

        title = tsv_file.readline().strip().split()
        title.append("gene_type")
        csv_writer.writerow(title)

        for line in tsv_file.readlines():
            info = line.split()
            g_name = info[1]

            if float(info[-1]) > 0:
                if g_name in name_type:
                    g_type = name_type[g_name]
                else:
                    g_name = g_name.lower()
                    if "trn" in g_name:
                        g_type = "tRNA"
                    elif "linc" in g_name:
                        g_type = "lncRNA"
                    elif "rp11" in g_name:
                        g_type = "protein_coding"
                    elif "snor" in g_name:
                        g_type = "snoRNA"
                    else:
                        g_type = crawler(g_name)
                        # g_type = "unknown_type"

                info.append(g_type)

                if g_type not in gene_type_count:
                    gene_type_count[g_type] = 1
                else:
                    gene_type_count[g_type] += 1

            csv_writer.writerow(info)
            my_genes.append(info)

    print(gene_type_count)
    return my_genes, gene_type_count


def mapping_tpm():
    name_type = {}
    for each in load_mapped_genes()[1:]:
        if float(each[-2]) > 0:
            g_name = each[0]
            g_type = each[-1]
            name_type[g_name] = g_type
        else:
            break

    tsv_file = open("./mapping_data/gene_tpm.tsv")
    my_genes = []

    with open("gene_tpm_with_type.csv", "w", newline="") as result_file:
        csv_writer = csv.writer(result_file, dialect="excel")

        title = tsv_file.readline().strip().split()
        title.append("gene_type")
        csv_writer.writerow(title)

        for line in tsv_file.readlines():
            info = line.split()
            g_name = info[1]

            if float(info[-1]) > 0:
                if g_name in name_type:
                    g_type = name_type[g_name]
                else:
                    g_name = g_name.lower()
                    if "trn" in g_name:
                        g_type = "tRNA"
                    elif "linc" in g_name:
                        g_type = "lncRNA"
                    elif "rp11" in g_name:
                        g_type = "protein_coding"
                    elif "snor" in g_name:
                        g_type = "snoRNA"
                    else:
                        g_type = crawler(g_name)
                        # g_type = "unknown_type"

                info.append(g_type)

            csv_writer.writerow(info)
            my_genes.append(info)
    print(len(my_genes))


def plot_gene_count(mode="sum"):
    """
        Available mode: sum, 1, 2, 3
    """
    my_genes = load_mapped_genes()
    gene_type_count = {}
    for each in my_genes[1:]:
        if mode == "sum":
            reads = float(each[-2])
        elif mode == "1":
            reads = float(each[-5])
        elif mode == "2":
            reads = float(each[-4])
        elif mode == "3":
            reads = float(each[-3])
        else:
            Exception("Unknown mode given")
            return
        if reads > 0:
            g_type = each[-1]
            if g_type not in gene_type_count:
                gene_type_count[g_type] = 1
            else:
                gene_type_count[g_type] += 1
        else:
            break

    to_plot = list(map(lambda x: [x, gene_type_count[x]], list(gene_type_count)))
    to_plot.sort(key=lambda x: x[1], reverse=True)

    plot_num = 20
    if mode == "sum":
        plt.title("Expressed gene per type for all samples (extremely low reads are not listed)")
    else:
        plt.title("gene counts per type for sample %s (extremely low reads are not listed)" % mode)
    plt.ylabel("gene types")
    plt.xlabel("gene counts in 3 samples")
    # plt.bar(x=[x[0] for x in to_plot[:plot_num]], height=[x[1] for x in to_plot[:plot_num]])
    plt.barh(y=[x[0] for x in to_plot[:plot_num]], width=[x[1] for x in to_plot[:plot_num]])
    plt.show()
    plt.pause(0)

    # my_genes.sort(key=lambda x: x[-2])
    print(len(my_genes))


def plot_reads_count(mode="sum"):
    """
        Available mode: sum, 1, 2, 3
    """
    my_genes = load_mapped_genes()

    gene_type_count = {}
    for each in my_genes[1:]:
        if mode == "sum":
            reads = float(each[-2])
        elif mode == "1":
            reads = float(each[-5])
        elif mode == "2":
            reads = float(each[-4])
        elif mode == "3":
            reads = float(each[-3])
        else:
            Exception("Unknown mode given")
            return
        if reads > 0:
            g_type = each[-1]
            if g_type not in gene_type_count:
                gene_type_count[g_type] = reads
            else:
                gene_type_count[g_type] += reads
        else:
            break

    to_plot = list(map(lambda x: [x, gene_type_count[x]], list(gene_type_count)))
    to_plot.sort(key=lambda x: x[1], reverse=True)

    plot_num = 20
    if mode == "sum":
        plt.title("Expressed gene Reads per type for all samples (extremely low reads are not listed)")
    else:
        plt.title("gene counts per type for sample %s (extremely low reads are not listed)" % mode)
    plt.ylabel("gene types")
    plt.xlabel("gene counts in 3 samples")
    # plt.bar(x=[x[0] for x in to_plot[:plot_num]], height=[x[1] for x in to_plot[:plot_num]])
    plt.barh(y=[x[0] for x in to_plot[:plot_num]], width=[x[1] for x in to_plot[:plot_num]])
    plt.show()
    plt.pause(0)

    # my_genes.sort(key=lambda x: x[-2])
    print(len(my_genes))


if __name__ == '__main__':
    plot_gene_count(mode="sum")
    plot_reads_count()
    # mapping_tpm()
