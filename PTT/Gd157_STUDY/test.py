import matplotlib.pyplot as plt
import numpy as np

def plot_piecewise_constant(boundaries, values):
    """
    Plot a piecewise constant distribution.
    
    Parameters:
    boundaries (list or array): List of n+1 boundary values defining the intervals.
    values (list or array): List of n values for each interval.
    """
    if len(boundaries) != len(values) + 1:
        raise ValueError("The number of boundaries must be one more than the number of values.")
    
    # Prepare data for step plot
    x = []
    y = []
    
    for i in range(len(values)):
        x.extend([boundaries[i], boundaries[i+1]])
        y.extend([values[i], values[i]])
    
    # Plot the piecewise constant function
    plt.step(x, y, where='post', color='b', linewidth=2)
    print(f"x = {x}")
    print(f"y = {y}")
    plt.xlabel('x')
    plt.ylabel('Value')
    plt.title('Piecewise Constant Distribution')
    plt.grid(True)
    plt.show()

# Example usage
boundaries = [0, 1, 3, 5]    # n+1 boundary values for n=3 intervals
values = [2, 4, 1]           # values for each interval
plot_piecewise_constant(boundaries, values)