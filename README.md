# Tutoriales de inicio para OpenSeesPy / OpenSeesPy Starting tutorials

## Descripción 
Este respositorio está dedicado a tutoriales de inicio con el OpenSeesPy.
Los tutoriales están pensados para pórticos de concreto reforzado modelados con elementos de fibras.
Bajo ninguna circunstancia estos tutoriales reemplazan un curso de análisis avanzado, donde se exponen fundamentos de modelación de manera clara y detallada. 
El propósito de estos tutoriales es ayudar a quienes se están iniciando en el mundo de OpenSeesPy y los tutoriales están basados en modelos de concreto reforzado. 
Estos tutoriales se basan en los ejemplos disponibles en el wiki de OpenSees: https://opensees.berkeley.edu/wiki/index.php?title=OpenSees_Examples_Manual_--_Structural_Models_%26_Anlyses. Aunque estos últimos están escritos en el intérprete de OpenSees para TCL, se recomienda visitarlos para ver muchas más aplicaciones de OpenSees que beneficiarán su aprendizaje.

## Requisitos de instalación / Requirements
Los tutoriales están escritos en su gran mayoría en Jupyter Notebooks. El primer tutorial indica qué deberás tener instalado para los tutoriales y proporciona los comandos para hacerlo. 
Estos tutoriales parten de la base de tener instalados OpenSeesPy, opsvis, vfo y opseestools. 

Esta última es una librería de rutinas de análisis que he desarrollado conjuntamente con mis estudiantes. La librería se puede instalar directamente desde PyPi: https://pypi.org/project/opseestools/ y el repositorio está disponible en: https://github.com/odarroyo/opseestools

## Descripción de los tutoriales / Tutorial description
Los tutoriales 2 al 7 son una progresión que inicia con un modelo elástico de un elemento en voladizo y que culmina en el análisis de un pórtico de concreto modelado con elementos de fibra, sometido a un análisis dinámico. Estos tutoriales están basados en el ejemplo disponible en el wiki de OpenSees: https://opensees.berkeley.edu/wiki/index.php?title=OpenSees_Example_5._2D_Frame,_3-story_3-bay,_Reinforced-Concrete_Section_%26_Steel_W-Section

Los tutoriales 8 y 9 se separan de los ejemplos del wiki de OpenSees, ya que están mas enfocados en la integración de OpenSeesPy con Python:
El tutorial 8 muestra cómo tomar ventaja de Python para aplicar un análisis dinámico incremental con un registro sísmico.
El tutorial 9 es un ejemplo de la integración de OpenSeesPy con la librería de multiprocesamiento de Python, aplicando un IDA de múltiples registros a un modelo no lineal con dos muros de concreto reforzado. Este puede ser modificado si se desea aplicarlo a cualquier otro tipo de modelo.

## Description
This repository is dedicated to introductory tutorials on OpenSeesPy.
The tutorials are designed for reinforced concrete frames modeled with fiber elements.
Under no circumstances do these tutorials replace an advanced analysis course, where modeling fundamentals are presented in a clear and detailed manner.
The goal of this tutorials is helping those who are starting with OpenSeesPy. They are based on reinforced concrete models.
The tutorials are based on the examples available in the OpenSees wiki: https://opensees.berkeley.edu/wiki/index.php?title=OpenSees_Examples_Manual_--_Structural_Models_%26_Anlyses. Although these are coded in the TCL interpreter, it is recommended to visit the website to learn several other applications of OpenSees that will be beneficial for your learning process.

## Requirements

The first tutorial indicates what you will need to have installed for the tutorials and provides the commands to do so.
These tutorials assume that OpenSeesPy, opsvis, vfo, and opseestools are installed.
The latter is a library of analysis routines that I have developed together with my students. The library can be installed directly from PyPi: https://pypi.org/project/opseestools/ and the repository is available at: https://github.com/odarroyo/opseestools

##  Tutorial description

Tutorials 2 to 7 are a progression starting with an elastic model of a cantilever element and culminating in the analysis of a concrete  frame modeled with fiber elements, subjected to a dynamic analysis. These tutorials are based on the example available on the OpenSees wiki: https://opensees.berkeley.edu/wiki/index.php?title=OpenSees_Example_5._2D_Frame,_3-story_3-bay,_Reinforced-Concrete_Section_%26_Steel_W-Section.

Tutorials 8 and 9 are separate from the OpenSees wiki examples, as they are more focused on the integration of OpenSeesPy with Python:
Tutorial 8 shows how to take advantage of Python to apply incremental dynamic analysis with a seismic log.
Tutorial 9 is an example of integrating OpenSeesPy with the Python multiprocessing library by applying a multi-log IDA to a nonlinear model with two reinforced concrete walls. This can be modified if you wish to apply it to any other type of model.
