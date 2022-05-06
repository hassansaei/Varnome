#!/usr/bin/env python
# coding: utf-8

# In[1]:


from setuptools import setup


# In[ ]:


setup(
    name='Singlepy',
    version='1.0.0',
    description='A tool to indentify and classify variants from WGS and WES data',
    author='Hassan Saei',
    author_email='hassan.saei@inserm.fr',
    license='AAAAAA',
    url='https://github.com/saei2021/Singlepy',
    packages=['singlepy'],
    package_dir={'singlepy': 'singlepy'},
    install_requires=['pandas', 'numpy', 'matplotlib', 'scikit-learn'],
)

