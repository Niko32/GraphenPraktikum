import pandas as pd 

# Pareto Optimierung

ratings = pd.DataFrame([
    ["Salted Caramel Muffin", 3, 2, 3, 3, 3], 
    ["Salted Caramel Kuchen", 4, 4, 4, 2, 3],
    ["Apfel Zimt Muffin", 3, 3, 2, 4, 3]
], columns=["Name", "Geschmack", "Mundgefühl", "Saftigkeit", "Preis", "Sättigung"])

print(ratings)