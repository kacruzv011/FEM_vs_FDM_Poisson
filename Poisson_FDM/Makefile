# ==========================================================
# Makefile para el Solucionador de Poisson con FDM
# Compila todas las versiones (serial y paralelas)
# ==========================================================

# Compilador
CXX = g++

# Banderas de compilación
CXXFLAGS = -O3 -Wall -std=c++17
OMPFLAGS = -fopenmp

# Directorios
SRCDIR = src
BINDIR = bin

# Lista de todos los programas a compilar
PROGS = \
    poisson_serial \
    poisson_parallel \
    poisson_collapse \
    poisson_sections \
    poisson_schedule \
    poisson_sincronizacion

TARGETS = $(addprefix $(BINDIR)/, $(PROGS))

# --- REGLAS ---
.PHONY: all clean

all: $(TARGETS)

# Regla de compilación genérica para todos los ejecutables.
# Todos necesitan la bandera de OpenMP porque todos usarán omp_get_wtime()
# para medir el tiempo de forma consistente.
$(TARGETS): $(BINDIR)/% : $(SRCDIR)/%.cpp
	@mkdir -p $(BINDIR)
	@echo "Compilando (FDM): $< -> $@"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

# Regla de limpieza
clean:
	@echo "Limpiando directorio de compilación 'Taller_OpenMP_Poisson'..."
	@rm -rf $(BINDIR) data fig
	@rm -f *.csv *.dat *.png *.gp