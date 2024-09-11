import cirq
import numpy as np

class BlockOperation(cirq.Operation):
    """class to subsume multiple gates within a block"""
    def __init__(self, operations):
        self.operations = operations

    def with_qubits(self, *new_qubits):
        if len(new_qubits) != len(self.qubits):
            raise ValueError("Incorrect number of qubits provided.")
        new_operations = [
            op.with_qubits(*new_qubits[:len(op.qubits)])
            for op in self.operations
        ]
        return BlockOperation(new_operations)

    @property
    def qubits(self):
        qubits = []
        for op in self.operations:
            qubits.extend(op.qubits)
        return tuple(set(qubits))  # Ensure no duplicate qubits

    def _decompose_(self):
        # Return a list of Moments to ensure operations are non-overlapping.
        return [cirq.Moment([op]) for op in self.operations]

    def __str__(self):
        return f"BlockOperation({', '.join(str(op) for op in self.operations)})"

    def __repr__(self):
        return f"BlockOperation({repr(self.operations)})"
    

class BlockGate(cirq.Gate):
    """Class to subsume multiple gates within a block."""
    
    def __init__(self, operations):
        self.operations = operations
        # Ensure all operations are of type cirq.Operation
        if not all(isinstance(op, cirq.Operation) for op in operations):
            raise ValueError("All items in 'operations' must be of type cirq.Operation")

        # Ensure all qubits are LineQubits
        for op in operations:
            for qubit in op.qubits:
                if not isinstance(qubit, cirq.LineQubit):
                    raise ValueError(f"Qubit {qubit} is not a LineQubit. Only LineQubits are allowed for the joint cutting with qsimh.")


    def with_qubits(self, *new_qubits):
        if len(new_qubits) != len(self.qubits):
            raise ValueError("Incorrect number of qubits provided.")
        new_operations = [
            op.with_qubits(*new_qubits[:len(op.qubits)])
            for op in self.operations
        ]
        return BlockGate(new_operations)

    @property
    def qubits(self):
        qubits = []
        for op in self.operations:
            qubits.extend(op.qubits)
        return tuple(set(qubits))  # Ensure no duplicate qubits

    def num_qubits(self):
        return len(set(qubit for op in self.operations for qubit in op.qubits))

    def _decompose_(self):
        # Return a list of Moments to ensure operations are non-overlapping.
        return [cirq.Moment([op]) for op in self.operations]

    def __str__(self):
        return f"BlockGate({', '.join(str(op) for op in self.operations)})"

    def __repr__(self):
        return f"BlockGate({repr(self.operations)})"

