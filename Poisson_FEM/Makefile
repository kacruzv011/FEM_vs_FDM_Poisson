# ==========================================================
# Makefile para el Solucionador de Poisson con FEM (Corregido)
# Compila todas las versiones.
# ==========================================================

# --- CONFIGURACIÓN DE COMPILACIÓN ---

# Compilador a utilizar
CXX = g++

# Banderas generales de optimización y advertencias
CXXFLAGS = -O3 -Wall -std=c++17

# Banderas específicas para OpenMP
OMPFLAGS = -fopenmp

# Ruta a la librería Eigen
EIGEN_INCLUDE = -I$(HOME)/libs/eigen

# --- CONFIGURACIÓN DE ARCHIVOS Y DIRECTORIOS ---

# Directorio donde se encuentran los archivos fuente (.cpp)
SRCDIR = src

# Directorio donde se guardarán los ejecutables compilados
BINDIR = bin

# Lista de todos los programas a compilar
PROGS = \
    poisson_serial \
    poisson_parallel \
    poisson_collapse \
    poisson_sections \
    poisson_schedule \
    poisson_sincronizacion

# Genera la lista completa de rutas a los ejecutables
TARGETS = $(addprefix $(BINDIR)/, $(PROGS))

# --- REGLAS DE MAKE ---

.PHONY: all clean

all: $(TARGETS)

# --- REGLAS DE COMPILACIÓN ---
# He simplificado las reglas para que sean más fáciles de leer.

# Regla para todos los ejecutables (ahora todos usan OpenMP para timing o paralelismo)
$(TARGETS): $(BINDIR)/% : $(SRCDIR)/%.cpp
	@mkdir -p $(BINDIR)
	@echo "Compilando: $< -> $@"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(EIGEN_INCLUDE) -o $@ $<

# --- REGLA DE LIMPIEZA ---

clean:
	@echo "Limpiando directorio de compilación..."
	@rm -rf $(BINDIR)