// Il programma del corso, nell'ordine delle lezioni. "simbolo" è la sigla della casella nella tavola del programma.
// pronto: l'argomento ha la sua pagina completa (teoria, esercizi, quiz).

export const PARTI = [
  { id: 'I', titolo: 'Parte I', sottotitolo: 'Principi di chimica generale' },
  { id: 'II', titolo: 'Parte II', sottotitolo: 'Equilibri, struttura e chimica inorganica' },
];

export const PROGRAMMA = [
  { n: 1, id: 'materia', parte: 'I', ore: 2, simbolo: 'Ma', titolo: 'Materia – Proprietà e Misura', pronto: true,
    fonte: 'Petrucci, <i>Chimica generale</i>, cap. 1 (Materia: proprietà e misura)' },
  { n: 2, id: 'atomi', parte: 'I', ore: 2, simbolo: 'At', titolo: 'Atomi e Teoria Atomica' },
  { n: 3, id: 'composti', parte: 'I', ore: 2, simbolo: 'Cp', titolo: 'Composti Chimici' },
  { n: 4, id: 'reazioni', parte: 'I', ore: 2, simbolo: 'Rz', titolo: 'Reazioni Chimiche' },
  { n: 5, id: 'soluzione-acquosa', parte: 'I', ore: 3, simbolo: 'Aq', titolo: 'Reazioni in Soluzione Acquosa' },
  { n: 6, id: 'gas', parte: 'I', ore: 2, simbolo: 'Gs', titolo: 'Gas Ideali e Reali' },
  { n: 7, id: 'primo-principio', parte: 'I', ore: 3, simbolo: 'Td', titolo: 'Primo Principio della Termodinamica con Applicazioni a Sistemi Chimici' },
  { n: 8, id: 'cinetica', parte: 'I', ore: 2, simbolo: 'Ci', titolo: 'Cinetica Chimica' },
  { n: 9, id: 'equilibrio', parte: 'I', ore: 2, simbolo: 'Eq', titolo: 'Equilibrio Chimico – Introduzione' },
  { n: 10, id: 'teoria-atomica-moderna', parte: 'I', ore: 3, simbolo: 'Qm', titolo: 'Teoria Atomica Moderna' },
  { n: 11, id: 'periodicita', parte: 'I', ore: 2, simbolo: 'Pe', titolo: 'Periodicità Chimica, Genesi degli Elementi e Radioattività' },
  { n: 12, id: 'legame', parte: 'I', ore: 6, simbolo: 'Lg', titolo: 'Legame Chimico – Modelli Molecolari e Strutturali' },
  { n: 13, id: 'liquidi-solidi', parte: 'I', ore: 3, simbolo: 'Ls', titolo: 'Liquidi, Solidi e Diagrammi di Stato – Cenni' },
  { n: 14, id: 'soluzioni', parte: 'I', ore: 2, simbolo: 'So', titolo: 'Soluzioni – Proprietà Generali' },
  { n: 15, id: 'acidi-basi', parte: 'II', ore: 4, simbolo: 'Ab', titolo: 'Equilibrio Chimico – Acidi e Basi' },
  { n: 16, id: 'solubilita', parte: 'II', ore: 2, simbolo: 'Ks', titolo: 'Equilibrio Chimico – Solubilità e Complessazione in Soluzione Acquosa' },
  { n: 17, id: 'gibbs', parte: 'II', ore: 3, simbolo: 'Gb', titolo: 'Entropia e Funzione di Gibbs' },
  { n: 18, id: 'elettrochimica', parte: 'II', ore: 3, simbolo: 'El', titolo: 'Elettrochimica' },
  { n: 19, id: 'simmetria', parte: 'II', ore: 6, simbolo: 'Sy', titolo: 'Simmetria Molecolare ed Elementi della Teoria dei Gruppi – Introduzione' },
  { n: 20, id: 'mo-lcao', parte: 'II', ore: 6, simbolo: 'Mo', titolo: 'Metodo MO-LCAO – Struttura e Reattività' },
  { n: 21, id: 'stato-solido', parte: 'II', ore: 2, simbolo: 'Ss', titolo: 'Stato Solido' },
  { n: 22, id: 'donatore-accettore', parte: 'II', ore: 2, simbolo: 'Da', titolo: 'Sistemi Acido-Base e Sistemi Donatore-Accettore' },
  { n: 23, id: 'gruppi-principali', parte: 'II', ore: 16, simbolo: 'Gp', titolo: 'Elementi dei Gruppi Principali' },
];

export const argomento = (id) => PROGRAMMA.find(a => a.id === id);
