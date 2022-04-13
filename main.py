import pandas as pd
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import requests
import gzip
import shutil
import os
import time
from zipfile import ZipFile

# working directory
os.chdir("/Users/USER/Desktop/NGS/")
print(os.getcwd())


def CADD_hom_input(csv_file='query.output.genome_summary.csv'):
    df = pd.read_csv(csv_file, low_memory=False)
    df = df.iloc[:, 0:135]
    dff_ = filter_on_zygosity_hom(df)
    df_ = filter_on_population(dff_)
    df_ = df_.iloc[:, 129:134]
    df_ = df_.rename(columns={"Otherinfo.3": "CHROM", "Otherinfo.4": "POS", "Otherinfo.5": "ID", "Otherinfo.6": "REF",
                              "Otherinfo.7": "ALT"})
    df_['ID'] = '.'
    write_vcf(df_, 'CADD_hom.vcf')
    return df, df_


def CADD_het_input(csv_file='query.output.genome_summary.csv'):
    df = pd.read_csv(csv_file, low_memory=False)
    df = df.iloc[:, 0:135]
    df2_ = filter_on_zygosity_het(df)
    df1_ = filter_on_population(df2_)
    df1_ = df1_.iloc[:, 129:134]
    df1_ = df1_.rename(columns={"Otherinfo.3": "CHROM", "Otherinfo.4": "POS", "Otherinfo.5": "ID", "Otherinfo.6": "REF",
                                "Otherinfo.7": "ALT"})
    df1_['ID'] = '.'
    write_vcf(df1_, 'CADD_het.vcf')
    return df, df1_


def ClinVar_results_hom(csv_file='query.output.genome_summary.csv'):
    df3 = pd.read_csv(csv_file, low_memory=False)
    df3 = df3.iloc[:, 0:135]
    df3 = filter_on_zygosity_hom(df3)
    df3 = filter_on_population(df3)
    df3_ = df3[df3['ClinVar_SIG'].str.contains("Pathogenic", "Likely pathogenic", na=False)]
    df3_.to_csv('ClinVar_Result_hom.csv', index=False)
    return df3_


def ClinVar_results_het(csv_file='query.output.genome_summary.csv'):
    df4 = pd.read_csv(csv_file, low_memory=False)
    df4 = df4.iloc[:, 0:135]
    df4 = filter_on_zygosity_het(df4)
    df4 = filter_on_population(df4)
    df4_ = df4[df4['ClinVar_SIG'].str.contains("Pathogenic", "Likely pathogenic", na=False)]
    df4_.to_csv('ClinVar_Result_het.csv', index=False)
    return df4_


def Functional_results_hom(csv_file='query.output.genome_summary.csv'):
    df = pd.read_csv(csv_file, low_memory=False)
    df5 = df.iloc[:, 0:135]
    df5 = filter_on_zygosity_hom(df5)
    df5 = filter_on_population(df5)
    df5_ = df5[df5['ExonicFunc.refGene'].isin(['frameshift deletion', 'frameshift insertion',
                                               'stopgain', 'startloss', 'stoploss'])]
    df5_.to_csv('Functional_Result_hom.csv', index=False)
    return df5_


def Functional_results_het(csv_file='query.output.genome_summary.csv'):
    df = pd.read_csv(csv_file, low_memory=False)
    df6 = df.iloc[:, 0:135]
    df6 = filter_on_zygosity_het(df6)
    df6 = filter_on_population(df6)
    df6_ = df6[df6['ExonicFunc.refGene'].isin(['frameshift deletion', 'frameshift insertion',
                                               'stopgain', 'startloss', 'stoploss'])]
    df6_.to_csv('Functional_Result_het.csv', index=False)
    return df6_


def filter_on_population(
        df,
        th_1000G_ALL=0.05,
        th_1000G_AFR=0.05,
        th_1000G_AMR=0.05,
        th_1000G_EAS=0.05,
        th_1000G_EUR=0.05,
        th_1000G_SAS=0.05,
        th_ExAC_Freq=0.05,
        th_ExAC_AFR=0.05,
        th_ExAC_AMR=0.05,
        th_ExAC_EAS=0.05,
        th_ExAC_FIN=0.05,
        th_ExAC_NFE=0.05,
        th_ExAC_OTH=0.05,
        th_ExAC_SAS=0.05,
        th_ESP6500si_ALL=0.05,
        th_ESP6500si_AA=0.05,
        th_ESP6500si_EA=0.05,
        th_gnomAD_exome_ALL=0.05,
        th_gnomAD_exome_AFR=0.05,
        th_gnomAD_exome_AMR=0.05,
        th_gnomAD_exome_ASJ=0.05,
        th_gnomAD_exome_EAS=0.05,
        th_gnomAD_exome_FIN=0.05,
        th_gnomAD_exome_NFE=0.05,
        th_gnomAD_exome_OTH=0.05,
        th_gnomAD_exome_SAS=0.05,
        th_gnomAD_genome_ALL=0.05,
        th_gnomAD_genome_AFR=0.05,
        th_gnomAD_genome_AMR=0.05,
        th_gnomAD_genome_ASJ=0.05,
        th_gnomAD_genome_EAS=0.05,
        th_gnomAD_genome_FIN=0.05,
        th_gnomAD_genome_NFE=0.05,
        th_gnomAD_genome_OTH=0.05
):
    df1_ = df.copy()
    df2_ = df1_.apply(lambda x: x.replace('.', -1.))
    df_ = df2_.fillna(0)

    df_['1000G_ALL'] = pd.to_numeric(df_['1000G_ALL'])
    df_['1000G_AFR'] = pd.to_numeric(df_['1000G_AFR'])
    df_['1000G_AMR'] = pd.to_numeric(df_['1000G_AMR'])
    df_['1000G_EAS'] = pd.to_numeric(df_['1000G_EAS'])
    df_['1000G_EUR'] = pd.to_numeric(df_['1000G_EUR'])
    df_['1000G_SAS'] = pd.to_numeric(df_['1000G_SAS'])
    df_['ExAC_Freq'] = pd.to_numeric(df_['ExAC_Freq'])
    df_['ExAC_AFR'] = pd.to_numeric(df_['ExAC_AFR'])
    df_['ExAC_AMR'] = pd.to_numeric(df_['ExAC_AMR'])
    df_['ExAC_EAS'] = pd.to_numeric(df_['ExAC_EAS'])
    df_['ExAC_FIN'] = pd.to_numeric(df_['ExAC_FIN'])
    df_['ExAC_NFE'] = pd.to_numeric(df_['ExAC_NFE'])
    df_['ExAC_OTH'] = pd.to_numeric(df_['ExAC_OTH'])
    df_['ExAC_SAS'] = pd.to_numeric(df_['ExAC_SAS'])
    df_['ESP6500si_ALL'] = pd.to_numeric(df_['ESP6500si_ALL'])
    df_['ESP6500si_AA'] = pd.to_numeric(df_['ESP6500si_AA'])
    df_['ESP6500si_EA'] = pd.to_numeric(df_['ESP6500si_EA'])
    df_['gnomAD_exome_ALL'] = pd.to_numeric(df_['gnomAD_exome_ALL'])
    df_['gnomAD_exome_AFR'] = pd.to_numeric(df_['gnomAD_exome_AFR'])
    df_['gnomAD_exome_AMR'] = pd.to_numeric(df_['gnomAD_exome_AMR'])
    df_['gnomAD_exome_ASJ'] = pd.to_numeric(df_['gnomAD_exome_ASJ'])
    df_['gnomAD_exome_EAS'] = pd.to_numeric(df_['gnomAD_exome_EAS'])
    df_['gnomAD_exome_FIN'] = pd.to_numeric(df_['gnomAD_exome_FIN'])
    df_['gnomAD_exome_NFE'] = pd.to_numeric(df_['gnomAD_exome_NFE'])
    df_['gnomAD_exome_OTH'] = pd.to_numeric(df_['gnomAD_exome_OTH'])
    df_['gnomAD_exome_SAS'] = pd.to_numeric(df_['gnomAD_exome_SAS'])
    df_['gnomAD_genome_ALL'] = pd.to_numeric(df_['gnomAD_genome_ALL'])
    df_['gnomAD_genome_AFR'] = pd.to_numeric(df_['gnomAD_genome_AFR'])
    df_['gnomAD_genome_AMR'] = pd.to_numeric(df_['gnomAD_genome_AMR'])
    df_['gnomAD_genome_ASJ'] = pd.to_numeric(df_['gnomAD_genome_ASJ'])
    df_['gnomAD_genome_EAS'] = pd.to_numeric(df_['gnomAD_genome_EAS'])
    df_['gnomAD_genome_FIN'] = pd.to_numeric(df_['gnomAD_genome_FIN'])
    df_['gnomAD_genome_NFE'] = pd.to_numeric(df_['gnomAD_genome_NFE'])
    df_['gnomAD_genome_OTH'] = pd.to_numeric(df_['gnomAD_genome_OTH'])

    df_ = df_[df_['1000G_ALL'] <= th_1000G_ALL]
    df_ = df_[df_['1000G_AFR'] <= th_1000G_AFR]
    df_ = df_[df_['1000G_AMR'] <= th_1000G_AMR]
    df_ = df_[df_['1000G_EAS'] <= th_1000G_EAS]
    df_ = df_[df_['1000G_EUR'] <= th_1000G_EUR]
    df_ = df_[df_['1000G_SAS'] <= th_1000G_SAS]
    df_ = df_[df_['ExAC_Freq'] <= th_ExAC_Freq]
    df_ = df_[df_['ExAC_AFR'] <= th_ExAC_AFR]
    df_ = df_[df_['ExAC_AMR'] <= th_ExAC_AMR]
    df_ = df_[df_['ExAC_EAS'] <= th_ExAC_EAS]
    df_ = df_[df_['ExAC_FIN'] <= th_ExAC_FIN]
    df_ = df_[df_['ExAC_NFE'] <= th_ExAC_NFE]
    df_ = df_[df_['ExAC_OTH'] <= th_ExAC_OTH]
    df_ = df_[df_['ExAC_SAS'] <= th_ExAC_SAS]
    df_ = df_[df_['ESP6500si_ALL'] <= th_ESP6500si_ALL]
    df_ = df_[df_['ESP6500si_AA'] <= th_ESP6500si_AA]
    df_ = df_[df_['ESP6500si_EA'] <= th_ESP6500si_EA]
    df_ = df_[df_['gnomAD_exome_ALL'] <= th_gnomAD_exome_ALL]
    df_ = df_[df_['gnomAD_exome_AFR'] <= th_gnomAD_exome_AFR]
    df_ = df_[df_['gnomAD_exome_AMR'] <= th_gnomAD_exome_AMR]
    df_ = df_[df_['gnomAD_exome_ASJ'] <= th_gnomAD_exome_ASJ]
    df_ = df_[df_['gnomAD_exome_EAS'] <= th_gnomAD_exome_EAS]
    df_ = df_[df_['gnomAD_exome_FIN'] <= th_gnomAD_exome_FIN]
    df_ = df_[df_['gnomAD_exome_NFE'] <= th_gnomAD_exome_NFE]
    df_ = df_[df_['gnomAD_exome_OTH'] <= th_gnomAD_exome_OTH]
    df_ = df_[df_['gnomAD_exome_SAS'] <= th_gnomAD_exome_SAS]
    df_ = df_[df_['gnomAD_genome_ALL'] <= th_gnomAD_genome_ALL]
    df_ = df_[df_['gnomAD_genome_AFR'] <= th_gnomAD_genome_AFR]
    df_ = df_[df_['gnomAD_genome_AMR'] <= th_gnomAD_genome_AMR]
    df_ = df_[df_['gnomAD_genome_ASJ'] <= th_gnomAD_genome_ASJ]
    df_ = df_[df_['gnomAD_genome_EAS'] <= th_gnomAD_genome_EAS]
    df_ = df_[df_['gnomAD_genome_FIN'] <= th_gnomAD_genome_FIN]
    df_ = df_[df_['gnomAD_genome_NFE'] <= th_gnomAD_genome_NFE]
    df_ = df_[df_['gnomAD_genome_OTH'] <= th_gnomAD_genome_OTH]

    df_ = df_.apply(lambda x: x.replace(-1., '.'))

    return df_


def filter_on_zygosity_het(df, th_Otherinfo='het'):
    df_ = df.copy()
    df2_ = df_[df_['Otherinfo'] == th_Otherinfo]

    return df2_


def filter_on_zygosity_hom(df, th_Otherinfo='hom'):
    df_ = df.copy()
    dff_ = df_[df_['Otherinfo'] == th_Otherinfo]

    return dff_


def write_vcf(df, output_VCF):
    with open(output_VCF, 'w') as vcf:
        df.to_csv(output_VCF, sep='\t', index=False)
        vcf.close()


CADD_hom_input()
CADD_het_input()
ClinVar_results_hom()
ClinVar_results_het()
Functional_results_hom()
Functional_results_het()


def request_to_CADD_hom(name='CADD_hom.vcf'):
    cwd = os.getcwd().replace("\\", "/")
    os.chdir(cwd)
    driver = webdriver.Chrome(executable_path="/Users/USER/Desktop/NGS/chromedriver/chromedriver.exe")
    driver.get("https://cadd.gs.washington.edu/score")
    f = driver.find_element_by_name("file")
    f.send_keys(cwd + "/" + name)
    version = driver.find_element_by_name("version")
    version.send_keys("GRCh37-v1.4")
    inc = driver.find_element_by_name("inclAnno")
    inc.click()
    driver.find_element_by_css_selector('form').submit()

    WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.ID, "content"), "successfully"))
    href = driver.find_element_by_css_selector('#content p:nth-child(2) a').get_attribute("href")
    driver.close()
    if "check_avail" in href:
        href = href.replace("check_avail", "static/finished")
    file_name = os.path.basename(href)

    r = requests.get(href)
    while r.status_code == 404:
        r = requests.get(href)
        time.sleep(5)

    with open(file_name, "wb") as file:
        file.write(r.content)

    n_file_name = file_name[:-3]
    with gzip.open(file_name, 'rb') as f_in:
        with open(n_file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return n_file_name


def request_to_CADD_het(name='CADD_het.vcf'):
    cwd = os.getcwd().replace("\\", "/")
    os.chdir(cwd)
    driver = webdriver.Chrome(executable_path="/Users/USER/Desktop/NGS/chromedriver/chromedriver.exe")
    driver.get("https://cadd.gs.washington.edu/score")
    f = driver.find_element_by_name("file")
    f.send_keys(cwd + "/" + name)
    version = driver.find_element_by_name("version")
    version.send_keys("GRCh37-v1.4")
    inc = driver.find_element_by_name("inclAnno")
    inc.click()
    driver.find_element_by_css_selector('form').submit()

    WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.ID, "content"), "successfully"))
    href = driver.find_element_by_css_selector('#content p:nth-child(2) a').get_attribute("href")
    driver.close()
    if "check_avail" in href:
        href = href.replace("check_avail", "static/finished")
    file_name = os.path.basename(href)

    r = requests.get(href)
    while r.status_code == 404:
        r = requests.get(href)
        time.sleep(5)

    with open(file_name, "wb") as file:
        file.write(r.content)

    n_file_name = file_name[:-3]
    with gzip.open(file_name, 'rb') as f_in:
        with open(n_file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return n_file_name


request_to_CADD_hom()
request_to_CADD_het()


def filter_base1_on_PHRED_14(d1_, th_PHRED=14):
    return d1_[d1_['PHRED'] >= th_PHRED]


def filter_base2_on_PHRED_14(d2_, th_PHRED=14):
    return d2_[d2_['PHRED'] >= th_PHRED]


def prepare_hom_for_mutationtaster(name=input("Please input HOM-TSV File :")):
    with open(name, 'r') as f:
        data = f.read().splitlines(True)
        f.close()
    with open(name, 'w') as f:
        f.writelines(data[1:])
        f.close()
    d1_ = pd.read_csv(name, sep='\t', low_memory=False)
    d1_ = filter_base1_on_PHRED_14(d1_)
    d1_ = d1_.iloc[:, 0:4]
    d1_.insert(2, 'ID', '.')
    d1_['.1'] = '.'
    d1_['.2'] = '.'
    d1_['.3'] = '.'
    d1_['.4'] = '.'
    d1_['.5'] = '.'
    write_vcf(d1_, 'hom_mut_taster.vcf')
    return d1_


def prepare_het_for_mutationtaster(name=input("PLease input HET-TSV File: ")):
    with open(name, 'r') as f:
        data = f.read().splitlines(True)
        f.close()
    with open(name, 'w') as f:
        f.writelines(data[1:])
        f.close()
    d2_ = pd.read_csv(name, sep='\t', low_memory=False)
    d2_ = filter_base2_on_PHRED_14(d2_)
    d2_ = d2_.iloc[:, 0:4]
    d2_.insert(2, 'ID', '.')
    d2_['.1'] = '.'
    d2_['.2'] = '.'
    d2_['.3'] = '.'
    d2_['.4'] = '.'
    d2_['.5'] = '.'
    write_vcf(d2_, 'het_mut_taster.vcf')
    return d2_


prepare_hom_for_mutationtaster()
prepare_het_for_mutationtaster()


def request_to_muttaster_hom(name=input("Please INPUT HOM-MUT VCF FILE:")):
    cwd = os.getcwd().replace("\\", "/")
    driver = webdriver.Chrome(executable_path="/Users/USER/Desktop/NGS/chromedriver/chromedriver.exe")
    driver.get("http://www.mutationtaster.org/StartQueryEngine.html")
    f = driver.find_element_by_name("filename")
    f.send_keys(cwd + "/" + name)
    project_name = driver.find_element_by_name("name")
    project_name.send_keys("HOM")
    driver.find_element_by_name("homo").click()
    driver.find_element_by_name("close_variants").click()
    driver.find_element_by_name("tgp_filter_hetero").click()
    driver.find_element_by_name("Submit").click()
    WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.CLASS_NAME, "border2"),
                                                                      "Display / filter / export results"))
    driver.find_element_by_class_name("border2").find_element_by_css_selector \
        ("table tbody tr:nth-child(2) td table tbody tr:nth-child(2) td:nth-child(2) input").click()
    driver.switch_to.window(driver.window_handles[1])
    WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.CLASS_NAME, "red"),
                                                                      "due to the recursion, "
                                                                      "limit/offset do not have any effect !"))
    href = driver.find_element_by_partial_link_text("retrieve TSV file").get_attribute("href")
    driver.quit()

    file_name = os.path.basename(href)
    r = requests.get(href)

    with open(file_name, "wb") as file:
        file.write(r.content)

    with ZipFile(file_name, "r") as file:
        file.extractall()

    return "HOM.tsv"


def request_to_muttaster_het(name=input("Please INPUT HET-MUT VCF FILE:")):
    cwd = os.getcwd().replace("\\", "/")
    driver = webdriver.Chrome(executable_path="/Users/USER/Desktop/NGS/chromedriver/chromedriver.exe")
    driver.get("http://www.mutationtaster.org/StartQueryEngine.html")
    f = driver.find_element_by_name("filename")
    f.send_keys(cwd + "/" + name)
    project_name = driver.find_element_by_name("name")
    project_name.send_keys("HET")
    driver.find_element_by_name("homo").click()
    driver.find_element_by_name("close_variants").click()
    driver.find_element_by_name("tgp_filter_hetero").click()
    driver.find_element_by_name("Submit").click()
    WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.CLASS_NAME, "border2"),
                                                                      "Display / filter / export results"))
    driver.find_element_by_class_name("border2").find_element_by_css_selector \
        ("table tbody tr:nth-child(2) td table tbody tr:nth-child(2) td:nth-child(2) input").click()
    driver.switch_to.window(driver.window_handles[1])
    WebDriverWait(driver, 600).until(EC.text_to_be_present_in_element((By.CLASS_NAME, "red"),
                                                                      "due to the recursion, "
                                                                      "limit/offset do not have any effect !"))
    href = driver.find_element_by_partial_link_text("retrieve TSV file").get_attribute("href")
    driver.quit()

    file_name = os.path.basename(href)
    r = requests.get(href)

    with open(file_name, "wb") as file:
        file.write(r.content)

    with ZipFile(file_name, "r") as file:
        file.extractall()

    return "HET.tsv"


request_to_muttaster_het()
request_to_muttaster_hom()


def prepare_final_result_hom(tsv_name_hom=input("Please add file for hom final prep:")):
    final_df_hom = pd.read_csv(tsv_name_hom, sep='\t')
    final_df_hom = final_df_hom[final_df_hom['pred_index'].isin(['disease_causing_automatic', 'disease_causing'])]
    final_df_hom['genesymbol'] = final_df_hom['genesymbol'].map(lambda x: x.rstrip('<br>compound'))
    final_df_hom = final_df_hom.drop_duplicates(subset=['genesymbol'])
    final_df_hom.to_csv('final_hom.csv', index=False)
    return final_df_hom


def prepare_final_result_het(tsv_name_het=input("Please add file for HET final prep:")):
    final_df_het = pd.read_csv(tsv_name_het, sep='\t')
    final_df_het = final_df_het[final_df_het['pred_index'].isin(['disease_causing_automatic', 'disease_causing'])]
    final_df_het['genesymbol'] = final_df_het['genesymbol'].map(lambda x: x.rstrip('<br>compound'))
    final_df_het = final_df_het.drop_duplicates(subset=['genesymbol'])
    final_df_het.to_csv('final_het.csv', index=False)
    return final_df_het


prepare_final_result_het()
prepare_final_result_hom()

