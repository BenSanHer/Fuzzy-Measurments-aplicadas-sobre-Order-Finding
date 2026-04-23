# Dirección 2 — Curvas de Información Mutua bajo Mediciones Difusas

## Motivación

Se estudia el efecto de mediciones difusas (fuzzy measurements) sobre el algoritmo de Shor, modeladas como canales que introducen permutaciones en el orden de los qubits medidos.

El objetivo es analizar cómo estas permutaciones afectan la información contenida en la distribución de probabilidad resultante.

## Modelo

Se considera una distribución ideal \( P_0(x) \) y una distribución perturbada \( P_p(x) \), donde:

- Con probabilidad \( p \), se aplica una permutación sobre los qubits
- Con probabilidad \( 1 - p \), se mide correctamente

Esto define un canal de la forma:

\[
P_p(x) = (1 - p)\, P_0(x) + p\, P_{\text{perm}}(x)
\]

donde \( P_{\text{perm}} \) depende del tipo de permutación.

## Tipos de permutaciones estudiadas

1. **Swaps entre dos qubits**
2. **Permutaciones circulares**
3. **Intercambio de todos los pares de qubits**

## Observación principal

Se estudia la información mutua entre la distribución ideal y la distribución perturbada:

\[
I(p) = I(P_0, P_p)
\]

y se analiza su comportamiento como función de \( p \in [0,1] \).

## Resultados observados

### Caso 1 — Swaps entre dos qubits

- La curva \( I(p) \) presenta:
  - Máximos en \( p = 0 \) y \( p = 1 \)
  - Mínimo en \( p = 0.5 \)

Esto es consistente con la simetría del canal:
- En \( p = 0 \): no hay perturbación
- En \( p = 1 \): la permutación es determinista (no hay mezcla)
- En \( p = 0.5 \): máxima incertidumbre

### Caso 2 — Permutaciones circulares

- Se observan curvas con forma similar, pero **no centradas en \( p = 0.5 \)**
- Los mínimos aparecen desplazados, típicamente alrededor de:
  \[
  p \approx 0.5 \pm 0.2
  \]

Esto indica una **ruptura de simetría** respecto al caso de swaps.

### Caso 3 — Intercambio de pares

- (En proceso de caracterización)

## Problema central

¿Por qué distintas permutaciones inducen curvas de información mutua con comportamientos cualitativamente diferentes?

En particular:

- ¿Qué propiedad de la permutación determina la posición del mínimo?
- ¿Por qué algunas permutaciones generan curvas simétricas y otras no?

## Hipótesis inicial

La forma de la curva \( I(p) \) depende de:

- La estructura de la permutación sobre los índices de los qubits
- El grado en que la permutación preserva o destruye correlaciones en la distribución de Shor

## Preguntas abiertas

- ¿Existe una clasificación de permutaciones según el comportamiento de \( I(p) \)?
- ¿Se puede predecir la posición del mínimo a partir de la estructura de la permutación?
- ¿Estas curvas contienen información suficiente para inferir el canal aplicado?
- ¿Qué relación existe entre estas curvas y la degradación del algoritmo de Shor?

## Estado actual

- Observación empírica de patrones en simulaciones
- Identificación de diferencias claras entre tipos de permutación

## Siguiente paso

- Formalizar matemáticamente la relación entre permutaciones y forma de \( I(p) \)
- Explorar posibles parametrizaciones del canal basadas en estas curvas