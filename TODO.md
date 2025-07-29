# TODO, mis notas

## Prácticas de verano

Lo principal es trabajar en el congreso al principio y en el paper despues, pero la principal tarea debería ser pasar el QAOA a un método propio utilizando QAOAAnsatz de Qiskit, y calcular la energia de manera directa $\bra{\beta,\gamma}H_c\ket{\beta,\gamma}$. Como tengo el hamiltoniano con Pauli Strings, puedo sacar directamente la matriz $H_c$. Utilizando esto debería funcionar bien QAOA.

Después de esto, sería ver si esto funciona en el Qmio, y si no, ver como hacer, pero puede ser pasar el circuito a otro lenguaje. Con esto, mas lo del paper, debería ser suficiente.

Los pasos serían:

- [x] Presentar TFG
- [x] Finalizar la presentación del congreso
- [x] Ir al congreso
- [] Pruebas en el Qmio
- [] Paper
  - [] Comparacion todo
  - [x] Comparacion exacto
  - [x] Comparacion arboles
  - [ ] Conclusion
  - [ ] Intro
- [x] Probar a resolver esto cambiando la matriz a la matriz de distancia y resolviendo el max-cut default. Creo que esto tendría sentido sobre todo para el QAOA.