

# Who hasn't been counted

In <- readNamedRegionFromFile(file = "countries.xlsx", name="In")

All <- readNamedRegionFromFile(file = "countries.xlsx", name="All")

Out <- All$Total [ (All$Total %in% In$In) ]

write.csv(Out, "out.csv")