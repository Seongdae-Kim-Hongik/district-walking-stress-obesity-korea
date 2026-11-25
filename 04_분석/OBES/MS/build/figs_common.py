# Shared figure helpers: numbers from numbers_obes.json (OBES_NUMBERS, default next to this file), half-up rounding, covariate labels and groups.
import os, json
from decimal import Decimal, ROUND_HALF_UP
HERE = os.path.dirname(os.path.abspath(__file__))
N = json.load(open(os.environ.get("OBES_NUMBERS", os.path.join(HERE, "numbers_obes.json")), encoding="utf-8"))
M = N["main"]
def r3(x):
    s = str(Decimal(repr(float(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))
    return str(Decimal(repr(float(x))).quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP)) if s == "1.000" else s
def hu(x, d=0, sign=False):
    s = str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)) if d else str(Decimal(repr(float(x))).quantize(Decimal(1), rounding=ROUND_HALF_UP))
    if s in ("-0", "-0.0", "-0.00", "-0.000"): s = s[1:]
    return ("+" + s if sign and not s.startswith("-") else s)
W, ED, ST, DR, RE, SM = "walking_practice_rate", "tertiary_educated_share", "perceived_stress_rate", "high_risk_drinking_rate_drinkers", "ext__restaurant_per10k", "current_smoking_rate"
LIV, FAC, INS, VIL, MV, GR, CAR = "basic_livelihood_rate", "public_sports_facilities_per_100k", "ext__instant_food_per10k", "village_sports_facilities_per_100k", "moderate_pa_rate", "eloc__dt_1yl20421e__1___1", "cars_per_capita"
LABEL = {W: "Walking", MV: "Moderate-to-vigorous physical activity", FAC: "Public sports facilities", VIL: "Neighborhood exercise facilities",
         GR: "Green zone area, log(1 + x)", CAR: "Registered cars (log)", RE: "General restaurants (log)", INS: "Instant-food outlets (log)",
         DR: "High-risk drinking among drinkers", SM: "Current smoking", ST: "Perceived stress", ED: "Tertiary education", LIV: "Basic livelihood security recipients",
         "urban_population_share": "Urban-zone population", "clinics_per_100k": "Clinics", "sex_ratio": "Sex ratio", "fiscal_independence": "Fiscal independence", "pop80_share": "Population aged 80 years or older"}
GROUPS = [("Physical activity", [W, MV]), ("Built environment", [FAC, VIL, GR, CAR]), ("Food environment", [RE, INS]),
          ("Behavioral and social conditions", [DR, SM, ST, ED, LIV]),
          ("Controls", ["urban_population_share", "clinics_per_100k", "sex_ratio", "fiscal_independence", "pop80_share"])]
