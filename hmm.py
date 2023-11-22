from collections import Counter
import numpy as np

def calculate_hmm_parameters(alignments):
    num_sequences = len(alignments)
    sequence_length = len(alignments[0])
    match_threshold = num_sequences / 2

    # Definir estados y transiciones
    states = calculate_states(alignments, match_threshold)
    transitions = calculate_transitions(states)

    # Calcular probabilidades de transición
    total_transitions = sum(transitions.values())
    transitions = {key: (value + 1) / (total_transitions + len(transitions)) for key, value in transitions.items()}

    # Calcular probabilidades de emisión
    emission_probabilities = calculate_emission_probabilities(alignments, states)

    # Probabilidades de background equiprobables
    symbols = set(''.join(alignments).replace('-', ''))
    background_probabilities = {symbol: 1 / len(symbols) for symbol in symbols} # Same background

    return states, transitions, emission_probabilities, background_probabilities

def calculate_states(alignments, match_threshold):
    return ['M' if sum(seq[col] != '-' for seq in alignments) > match_threshold else 'I'
            for col in range(len(alignments[0]))]

def calculate_transitions(states):
    return Counter(zip(states, states[1:]))

def calculate_emission_probabilities(alignments, states):
    emission_probabilities = {'M': Counter(), 'I': Counter()}

    for i, state in enumerate(states):
        for seq in alignments:
            if seq[i] != '-':
                emission_probabilities[state][seq[i]] += 1

    for state in emission_probabilities:
        total = sum(emission_probabilities[state].values())
        emission_probabilities[state] = {symbol: (count + 1) / (total + len(emission_probabilities[state]))
                                         for symbol, count in emission_probabilities[state].items()}

    return emission_probabilities

alignments = [
    "VGA--HAGEY",
    "V----NVDEV",
    "VEA--DVAGH",
    "VKG------D",
    "VYS--TYETS",
    "FNA--NIPKH",
    "IAGADNGAGY"
]
states, transitions, emission_probabilities, background_probabilities = calculate_hmm_parameters(alignments)

# Imprimir resultados de manera más organizada
print("States:", states)
print("\nTransitions:")
for key, value in transitions.items():
    print(f"  {key}: {value}")
print("\nEmission Probabilities:")
for state, probs in emission_probabilities.items():
    print(f"  {state}: {probs}")
print("\nBackground Probabilities:")
for symbol, prob in background_probabilities.items():
    print(f"  {symbol}: {prob}")
