# Análisis Comparativo y Paralelización de Solucionadores para la Ecuación de Poisson

Este repositorio contiene el código fuente, los scripts de automatización y los resultados de un estudio comparativo entre el Método de Diferencias Finitas (FDM) y el Método de Elementos Finitos (FEM) para resolver la ecuación de Poisson en 2D. El proyecto explora diversas estrategias de paralelización utilizando OpenMP y automatiza todo el flujo de trabajo experimental mediante Jupyter Notebooks.

## 📝 Descripción del Proyecto

El objetivo principal de este trabajo es doble:
1.  **Implementar y validar** dos familias de solucionadores numéricos (FDM y FEM) para la ecuación de Poisson, un problema fundamental en física e ingeniería.
2.  **Analizar y comparar** el rendimiento y la escalabilidad de múltiples estrategias de paralelización (serial, paralelo estándar, `collapse`, `sections`, `schedule`, etc.) aplicadas a cada método.

El proyecto está estructurado en dos sub-proyectos autocontenidos:
*   `Poisson_FDM`: Contiene las implementaciones basadas en el **Método de Diferencias Finitas (FDM)**.
*   `Poisson_FEM`: Contiene las implementaciones basadas en el **Método de Elementos Finitos (FEM)**.

Todo el proceso, desde la compilación y ejecución de benchmarks hasta el análisis de datos y la visualización de resultados, está orquestado mediante Makefiles y Jupyter Notebooks para garantizar la **reproducibilidad** y la **eficiencia** del análisis.

## 📂 Estructura de Carpetas

```
/
├── Taller_OpenMP_Poisson/      # --- Proyecto de Diferencias Finitas (FDM) ---
│   ├── src/                    # Código fuente C++ para FDM
│   ├── bin/                    # Ejecutables compilados de FDM
│   ├── Makefile                # Makefile para compilar los programas FDM
│   └── analisis_fdm.ipynb      # Notebook para automatizar y analizar FDM
│
├── Poisson_FEM/                # --- Proyecto de Elementos Finitos (FEM) ---
│   ├── src/                    # Código fuente C++ para FEM
│   ├── bin/                    # Ejecutables compilados de FEM
│   ├── Makefile                # Makefile para compilar los programas FEM
│   └── analisis_fem.ipynb      # Notebook para automatizar y analizar FEM
│
├── Makefile                    # Makefile maestro para compilar todo
├── comparacion_final.ipynb     # Notebook para comparar los resultados de FDM vs. FEM
├── informe_final.pdf           # Documento LaTeX con el informe final del proyecto
└── README.md                   # Este archivo
```

## 🛠️ Requisitos y Dependencias

Para compilar y ejecutar este proyecto, necesitas un entorno basado en Linux (se recomienda WSL - Windows Subsystem for Linux) con las siguientes herramientas instaladas:

### Compilación (C++)
*   **Compilador C++:** `g++` (versión 7 o superior)
*   **Make:** `make`
*   **OpenMP:** Generalmente incluido con `g++`.
*   **Librería Eigen:** Una librería de C++ para álgebra lineal, requerida para el proyecto FEM.
    *   *Instrucciones de instalación:* El `Makefile` de FEM asume que Eigen se encuentra en `~/libs/eigen`. Puedes instalarlo con los siguientes comandos en tu terminal WSL:
        ```bash
        cd ~
        wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        tar -xvf eigen-3.4.0.tar.gz
        mkdir -p ~/libs
        mv eigen-3.4.0 ~/libs/eigen
        ```

### Análisis y Automatización (Python)
*   **Python 3:** (versión 3.8 o superior)
*   **Entorno Virtual (Recomendado):**
    ```bash
    # Desde la raíz de tu proyecto
    python3 -m venv venv
    source venv/bin/activate
    ```
*   **Librerías de Python:** Instalar con pip dentro del entorno virtual:
    ```bash
    pip install jupyterlab pandas matplotlib seaborn tqdm numpy
    ```

## 🚀 Cómo Empezar

Sigue estos pasos para compilar y ejecutar todo el proyecto.

### 1. Compilar los Programas

Con el `Makefile` maestro en la raíz del proyecto, puedes compilar ambos subproyectos con un solo comando:

```bash
# Desde la carpeta raíz del proyecto
make
```
Si prefieres compilarlos por separado:
```bash
cd Taller_OpenMP_Poisson && make && cd ..
cd Poisson_FEM && make && cd ..
```
Tras la compilación, las carpetas `bin/` dentro de cada proyecto contendrán todos los ejecutables.

### 2. Ejecutar los Benchmarks y Análisis

El análisis se realiza a través de los Jupyter Notebooks.

1.  **Activa tu entorno virtual de Python:**
    ```bash
    # Desde la raíz del proyecto
    source venv/bin/activate 
    ```
2.  **Inicia Jupyter Lab:**
    ```bash
    jupyter lab
    ```
3.  **Ejecuta los notebooks en el siguiente orden:**
    *   **`Taller_OpenMP_Poisson/analisis_fdm.ipynb`**: Ejecutará el benchmark para FDM, generando `resultados_fdm.csv` y mostrando sus gráficas.
    *   **`Poisson_FEM/analisis_fem.ipynb`**: Hará lo mismo para FEM, generando `resultados_fem.csv`.
    *   **`comparacion_final.ipynb`**: Cargará los dos archivos CSV generados y creará las gráficas comparativas finales.

## 📈 Resumen de Resultados y Conclusiones

*[Esta es una sección muy importante. Rellénala con tus hallazgos clave después de ejecutar los análisis. Aquí tienes un ejemplo de cómo podrías redactarla.]*

El análisis comparativo entre los métodos FDM y FEM reveló varios puntos clave:

*   **Rendimiento Serial:** Para mallas pequeñas, el método FDM demostró ser más rápido debido a su menor sobrecarga computacional. Sin embargo, en mallas de mayor tamaño, el método FEM, gracias al uso de solucionadores directos de álgebra lineal dispersa de la librería Eigen, fue significativamente más eficiente que el método iterativo de Jacobi de FDM.

*   **Escalabilidad Paralela:** El método FDM mostró una escalabilidad casi ideal. Su algoritmo de Jacobi es inherentemente paralelo, con mínima comunicación y sincronización entre hilos. Por el contrario, las implementaciones de FEM mostraron una buena escalabilidad hasta un número moderado de hilos (ej. 8), pero su rendimiento se estancó o degradó con un mayor número de núcleos. Esto se atribuye a la contención generada por la sección crítica utilizada para ensamblar la matriz global.

*   **Conclusión General:** No hay un único "mejor" método; la elección depende de los requisitos del problema.
    *   **FEM** es superior en rendimiento para problemas grandes y complejos donde la precisión y el manejo de geometrías son importantes, siempre que se utilicen librerías de álgebra lineal optimizadas.
    *   **FDM** es más simple de implementar y muestra una mejor escalabilidad paralela, lo que lo hace una opción viable para problemas en dominios regulares donde se puede aprovechar masivamente el paralelismo de datos.

Este proyecto subraya la importancia de una elección cuidadosa tanto del algoritmo numérico como de la estrategia de paralelización para lograr un rendimiento óptimo en la computación científica.
