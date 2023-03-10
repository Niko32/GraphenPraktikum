import pandas as pd 
from matplotlib import pyplot as plt

# Pareto Optimierung

ratings = pd.DataFrame([
    [3, 2, 3, 3, 3],
    [4, 4, 4, 2, 3],
    [3, 3, 2, 4, 3],
    [5, 5, 4, 4, 5],
    [5, 5, 5, 2, 2],
    [5, 4, 4, 3, 3]
], columns=["Geschmack", "Mundgefühl", "Saftigkeit", "Preis/Menge", "Sättigung"], 
    index=["Salted Caramel Muffin", "Salted Caramel Kuchen", "Apfel Zimt Muffin", "Mohnkuchen", "Himbeer-Kaesekuchen", "Carrot Cheese Cake"])

x = ratings.transpose().plot()
plt.yticks(range(6))
plt.savefig("output/plots/mensa.png")