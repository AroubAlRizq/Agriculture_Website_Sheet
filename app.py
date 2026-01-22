from flask import Flask, render_template, request, jsonify, send_file, url_for
import math
from flask import Flask, render_template, request, jsonify, send_file
import math
import pandas as pd
import io

app = Flask(__name__)

# --- CONFIGURATION: Define requirements for each formula ---
FORMULA_MAP = {
    # Evapotranspiration
    'hargreaves': ['T_mean', 'T_max', 'T_min', 'Ra'],
    'blaney_criddle': ['Ta', 'Kc'],
    'fao56': ['Rn', 'G', 'T_mean', 'u2', 'es_ea', 'Delta', 'Gamma'],
    'stephens_stewart': ['Ta', 'Rl'],
    'grassi': ['Ta', 'Rl'], 
    'linarce': ['Ta', 'Td', 'z', 'lat'],

    # GDD
    'gdd_arnold': ['TM', 'Tm', 'Tb'],
    'gdd_villa_nova': ['TM', 'Tm', 'Tb'],
    'gdd_ometto': ['TM', 'Tm', 'Tb', 'TB'],
    'gdd_snyder': ['TM', 'Tm', 'Tb', 'TB'],

    # Chill Units
    'chill_utah': ['T_current'],
    'chill_nc': ['T_current'],

    # Vegetation Indices
    'ndvi': ['NIR', 'Red'],
    'gndvi': ['NIR', 'Green'],
    'pri': ['R531', 'R570'],
    'ndre': ['R790', 'R720'],
    'ccci': ['NDRE', 'NDRE_min', 'NDRE_max'],
    'rvi': ['NIR', 'Red'],
    'evi': ['NIR', 'Red', 'Blue'],
    'evi2': ['NIR', 'Red'],
    'varigreen': ['Green', 'Red', 'Blue'],
    'vari700': ['R700', 'Red', 'Blue'],
    'tvi': ['R750', 'R550', 'R670'],
    'mtvi1': ['R800', 'R550', 'R670'],
    'mtvi2': ['R800', 'R550', 'R670'],
    'mtci': ['R753', 'R708', 'R681'],
    'car': ['R700', 'R670', 'R550'],
    'cari': ['R700', 'R670', 'R550'],
    'tcari': ['R700', 'R670', 'R550'],
    'mcari': ['R700', 'R670', 'R550'],
    'mcari1': ['R800', 'R670', 'R550'],
    'mcari2': ['R800', 'R670', 'R550'],
    'wdvi': ['NIR', 'Red', 'a'],
    'pvi': ['NIR', 'Red', 'a', 'b'],
    'savi': ['NIR', 'Red', 'L'],
    'tsavi': ['R800', 'R670', 'a', 'b'],
    'osavi': ['NIR', 'Red'],
    'msavi': ['NIR', 'Red', 'L'],
    'msavi2': ['NIR', 'Red'],
    'sarvi': ['R800', 'Red', 'Blue'],
}

# --- CORE CALCULATION ENGINE ---
def compute_value(formula, inputs):
    result = 0
    unit = ""

    # Helper: Safe Float Conversion with Case Insensitivity
    # This works for both Pandas Series (Batch) and Dictionary (Manual)
    def val(key):
        try:
            # Direct lookup first (Fast for Manual)
            if key in inputs:
                return float(inputs[key])
            # Case insensitive lookup (Robust for Batch Excel headers)
            for k in inputs.keys():
                if k.strip().lower() == key.lower():
                    return float(inputs[k])
            return 0.0 
        except:
            return 0.0


    try:
        # --- Helper: Safe Float Conversion ---
        def val(key):
            return float(inputs.get(key, 0))

        # ==========================================
        # 1. [cite_start]EVAPOTRANSPIRATION (ET) [cite: 2, 5]
        # ==========================================
        if formula == "hargreaves":
            # [cite_start]ET = 0.0023 * (Tmean + 17.8) * sqrt(Tmax-Tmin) * Ra [cite: 2]
            t_mean, t_max, t_min, ra = val('T_mean'), val('T_max'), val('T_min'), val('Ra')
            result = 0.0023 * (t_mean + 17.8) * math.sqrt(max(0, t_max - t_min)) * ra
            unit = "mm/day"
        
        elif formula == "blaney_criddle":
            # [cite_start]PET = (0.0173*Ta - 0.314) * Kc * Ta [cite: 2] (Simplified adaptation)
            ta, kc = val('Ta'), val('Kc')
            result = (0.0173 * ta - 0.314) * kc * ta
            unit = "mm/day"

        elif formula == "fao56":
            # [cite_start]Penman-Monteith [cite: 2]
            # T_o = (0.408 Delta (Rn-G) + Gamma (900/(T+273)) u2 (es-ea)) / (Delta + Gamma(1+0.34u2))
            rn, g, t, u2 = val('Rn'), val('G'), val('T_mean'), val('u2')
            es_ea, delta, gamma = val('es_ea'), val('Delta'), val('Gamma')
            
            num = (0.408 * delta * (rn - g)) + (gamma * (900 / (t + 273)) * u2 * es_ea)
            den = delta + (gamma * (1 + 0.34 * u2))
            result = num / den if den != 0 else 0
            unit = "mm/day"

        elif formula == "stephens_stewart":
            # [cite_start]PET = (0.0082 Ta - 0.19) * (Rl/1500) * 25.4 [cite: 5]
            ta, rl = val('Ta'), val('Rl')
            result = (0.0082 * ta - 0.19) * (rl / 1500.0) * 25.4
            unit = "mm/day"

        elif formula == "grassi":
            # [cite_start]PET = 0.537 * 0.000675 * Rl * (0.62 + 0.00559 * Ta) * 25.4 [cite: 5]
            # Assumed Cloud cover factor is part of Rl or coeff. strictly following text structure:
            ta, rl = val('Ta'), val('Rl')
            result = 0.537 * 0.000675 * rl * (0.62 + 0.00559 * ta) * 25.4
            unit = "mm/day"

        elif formula == "linarce":
            # [cite_start]PET = [700 Tm / (100-lat) + 15(Ta-Td)] / (80-Ta) [cite: 5]
            # Tm = Ta + 0.006 * z
            ta, td, z, lat = val('Ta'), val('Td'), val('z'), val('lat')
            tm = ta + 0.006 * z
            num = ((700 * tm) / (100 - lat)) + (15 * (ta - td))
            den = 80 - ta
            result = num / den if den != 0 else 0
            unit = "mm/day"

        # ==========================================
        # 2. [cite_start]GROWING DEGREE DAYS (GDD) [cite: 11, 13]
        # ==========================================
        elif formula == "gdd_arnold":
            tm_big, tm_small, tb = val('TM'), val('Tm'), val('Tb')
            result = ((tm_big + tm_small) / 2) - tb
            if result < 0: result = 0
            unit = "°C-days"

        elif formula == "gdd_villa_nova":
            #[cite_start]# [cite: 11]
            TM, Tm, Tb = val('TM'), val('Tm'), val('Tb')
            if Tb >= TM: # to check if base emp already less than the max temp of the day
                result = 0 # no growth because the highest temp was already too cold for the plants to grow
            elif Tb < Tm:
                result = ((Tm - Tb) + (TM - Tm)) / 2.0 
                # Note: This mathematically simplifies to (TM+Tm)/2 - Tb (Arnold), but uses Source 11 syntax.
            elif Tb>=Tm: 
                # Case: Tm < Tb < TM
                result = ((TM - Tb)**2) / (2 * (TM - Tm)) if (TM-Tm) != 0 else 0

            else: result = 0

            unit = "°C-days"

        elif formula == "gdd_ometto":
            #[cite_start]# [cite: 11]
            TM, Tm, Tb, TB = val('TM'), val('Tm'), val('Tb'), val('TB')
            
            
            if TM > TB and TB > Tm and Tm > Tb: # Case 4
                term1 = 2 * (TM - Tm) * (Tm - Tb)
                term2 = (TM - Tm)**2
                term3 = (TM - TB)**2
                den = 2 * (TM - Tm)
                result = (term1 + term2 - term3) / den
            
            
            elif TB > TM and TM > Tm and Tm > Tb: # Case 1
                result = ((TM - Tm) / 2.0) + (Tm - Tb)
                
            
            elif TM > TB and TB > Tb and Tb > Tm: # Case 5
                term1 = (TM - Tb)**2
                term2 = (TM - TB)**2
                den = 2 * (TM - Tm)
                result = 0.5 * ((term1 - term2) / (TM - Tm))

            
            elif TB > TM and TM > Tb and Tb > Tm: # Case 2
                 result = ((TM - Tb)**2) / (2 * (TM - Tm))

            elif TB > Tb and Tb > TM and TM > Tm: # Case 3
                 result = 0
            
            else:
                result = 0 # Default/Fallback
            unit = "°C-days"

        elif formula == "gdd_snyder":
            #[cite_start]# [cite: 11, 13] Sine Method
            TM, Tm, Tb, TB = val('TM'), val('Tm'), val('Tb'), val('TB')
            M = (TM + Tm) / 2.0
            W = (TM - Tm) / 2.0
            
            # Helper to calculate theta/phi
            def get_angle(thresh):
                val = (thresh - M) / W
                val = max(-1, min(1, val)) # Clamp for arcsin
                return math.asin(val)

            theta = get_angle(Tb)
            phi = get_angle(TB)
            pi = math.pi

            # Calculate GDD1 (Above Tb)
            if Tb <= Tm:
                gdd1 = M - Tb
            elif Tb >= TM:
                gdd1 = 0
            else:
                # Tm < Tb < TM
                gdd1 = ((1/pi) * ((M - Tb) * (pi/2 - theta) + (W * math.cos(theta))))

            # Calculate GDD2 (Above TB - penalty)
            if TB <= Tm:
                gdd2 = M - TB
            elif TB >= TM:
                gdd2 = 0
            else:
                gdd2 = ((1/pi) * ((M - TB) * (pi/2 - phi) + (W * math.cos(phi))))

            result = gdd1 - gdd2
            unit = "°C-days"

        # ==========================================
        # 3. [cite_start]CHILL UNITS [cite: 7, 9]
        # ==========================================
        elif formula == "chill_utah":
            # Lookup Table based on Source 9
            t = val('T_current')
            if t <= 1.4: result = 0
            elif 1.5 <= t <= 2.4: result = 0.5
            elif 2.5 <= t <= 9.1: result = 1.0
            elif 9.2 <= t <= 12.4: result = 0.5
            elif 12.5 <= t <= 15.9: result = 0
            elif 16.0 <= t <= 18.0: result = -0.5
            elif t > 18.0: result = -1.0
            unit = "Chill Units (Hourly)"

        elif formula == "chill_nc":
            # North Carolina Model (Source 9 Lower Table)
            t = val('T_current')
            if t <= 1.5: result = 0
            elif 1.6 <= t <= 7.1: result = 0.5
            elif 7.2 <= t <= 12.9: result = 1.0
            elif 13.0 <= t <= 16.4: result = 0.5
            elif 16.5 <= t <= 19.0: result = 0 # Review with Dr. Alaa
            elif 19.1 <= t <= 20.6: result = -0.5
            elif 20.7 <= t <= 22.0: result = -1.0
            elif 22.1 <= t <= 23.2: result = -1.5
            else: result = -2.0
            unit = "Chill Units (Hourly)"

        # ==========================================
        # 4. [cite_start]VEGETATION INDICES [cite: 16, 18, 19]
        # ==========================================
        elif formula == "ndvi":
            nir, red = val('NIR'), val('Red')
            result = (nir - red) / (nir + red) if (nir+red)!=0 else 0
            unit = "Index"

        elif formula == "gndvi":
            nir, green = val('NIR'), val('Green')
            result = (nir - green) / (nir + green) if (nir+green)!=0 else 0
            unit = "Index"

        elif formula == "pri":
            r531, r570 = val('R531'), val('R570')
            result = (r531 - r570) / (r531 + r570) if (r531+r570)!=0 else 0
            unit = "Index"

        elif formula == "ndre":
            # (R790 - R720) [cite_start]/ (R790 + R720) [cite: 18]
            nir, red_edge = val('R790'), val('R720')
            result = (nir - red_edge) / (nir + red_edge) if (nir+red_edge)!=0 else 0
            unit = "Index"

        elif formula == "ccci":
            # (NDRE - NDRE_min) [cite_start]/ (NDRE_max - NDRE_min) [cite: 18]
            ndre, n_min, n_max = val('NDRE'), val('NDRE_min'), val('NDRE_max')
            result = (ndre - n_min) / (n_max - n_min) if (n_max-n_min)!=0 else 0
            unit = "Index"

        elif formula == "rvi":
            nir, red = val('NIR'), val('Red')
            result = nir / red if red!=0 else 0
            unit = "Ratio"

        elif formula == "evi":
            # [cite_start]2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1) [cite: 18]
            nir, red, blue = val('NIR'), val('Red'), val('Blue')
            den = nir + (6 * red) - (7.5 * blue) + 1
            result = 2.5 * ((nir - red) / den) if den!=0 else 0
            unit = "Index"

        elif formula == "evi2":
            # [cite_start]2.5 * (NIR - Red) / (NIR + 2.4*Red + 1) [cite: 18]
            nir, red = val('NIR'), val('Red')
            den = nir + (2.4 * red) + 1
            result = 2.5 * ((nir - red) / den) if den!=0 else 0
            unit = "Index"
        
        elif formula == "varigreen":
            # (Green - Red) [cite_start]/ (Green + Red - Blue) [cite: 18]
            green, red, blue = val('Green'), val('Red'), val('Blue')
            den = green + red - blue
            result = (green - red) / den if den!=0 else 0
            unit = "Index"
        
        elif formula == "vari700":
             # (R700 - 1.7*Red + 0.7*Blue) [cite_start]/ (R700 + 2.3*Red - 1.3*Blue) [cite: 18]
            r700, red, blue = val('R700'), val('Red'), val('Blue')
            num = r700 - (1.7 * red) + (0.7 * blue)
            den = r700 + (2.3 * red) - (1.3 * blue)
            result = num / den if den!=0 else 0
            unit = "Index"

        elif formula == "tvi":
            # [cite_start]0.5 * [120*(R750-R550) - 200*(R670-R550)] [cite: 18]
            r750, r550, r670 = val('R750'), val('R550'), val('R670')
            result = 0.5 * (120 * (r750 - r550) - 200 * (r670 - r550))
            unit = "Index"

        elif formula == "mtvi1":
            # [cite_start]1.2 * [1.2*(R800-R550) - 2.5*(R670-R550)] [cite: 18]
            r800, r550, r670 = val('R800'), val('R550'), val('R670')
            result = 1.2 * (1.2 * (r800 - r550) - 2.5 * (r670 - r550))
            unit = "Index"
            
        elif formula == "mtvi2":
            #[cite_start]# [cite: 18] Complex formula involving MTVI1 logic inside numerator
            r800, r550, r670 = val('R800'), val('R550'), val('R670')
            num = 1.5 * (1.2 * (r800 - r550) - 2.5 * (r670 - r550))
            den = math.sqrt((2 * r800 + 1)**2 - (6 * r800 - 5 * math.sqrt(r670)) - 0.5)
            result = num / den if den!=0 else 0
            unit = "Index"

        elif formula == "mtci":
            # (R753.75 - R708.75) [cite_start]/ (R708.75 - R681.25) [cite: 18]
            r753, r708, r681 = val('R753'), val('R708'), val('R681')
            den = r708 - r681
            result = (r753 - r708) / den if den!=0 else 0
            unit = "Index"
        
        elif formula == "car":
            # CAR [Kim et al. [cite_start]1994] [cite: 18]
            r700, r670, r550 = val('R700'), val('R670'), val('R550')
            a = (r700 - r550) / 150.0 
            b = r550 - (a * 550)    
            term1 = abs(a * 670 + b + r670)
            term2 = math.sqrt(a**2 + 1)
            result = term1/term2
            unit = "Index"

        elif formula == "cari":
            # CARI [Kim et al. [cite_start]1994] [cite: 18]
            r700, r670, r550 = val('R700'), val('R670'), val('R550')
            a = (r700 - r550) / 150.0 
            b = r550 - (a * 550)    
            term1 = abs(a * 670 + b + r670)
            term2 = math.sqrt(a**2 + 1)
            term3 = term1/term2
            term4 = r700/r670
            result = term3*term4 if r670!=0 else 0
            unit = "Index"

        elif formula == "tcari":
             # TCARI [Haboudane et al. [cite_start]2002] [cite: 19]
            r700, r670, r550 = val('R700'), val('R670'), val('R550')
            term1 = r700 - r670
            term2 = 0.2 * (r700 - r550) * (r700 / r670) if r670!=0 else 0
            result = 3 * (term1 - term2)
            unit = "Index"

        elif formula == "mcari":
            # [cite_start][(R700-R670) - 0.2*(R700-R550)] * (R700/R670) [cite: 19]
            r700, r670, r550 = val('R700'), val('R670'), val('R550')
            result = ((r700 - r670) - 0.2 * (r700 - r550)) * (r700 / r670) if r670!=0 else 0
            unit = "Index"

        elif formula == "mcari1":
             # [cite_start]1.2 * [2.5*(R800-R670) - 1.3*(R800-R550)] [cite: 19]
            r800, r670, r550 = val('R800'), val('R670'), val('R550')
            result = 1.2 * (2.5 * (r800 - r670) - 1.3 * (r800 - r550))
            unit = "Index"

        elif formula == "mcari2":
            #[cite_start]# [cite: 19]
            r800, r670, r550 = val('R800'), val('R670'), val('R550')           
            num = 1.5 * (2.5 * (r800 - r670) - 1.3 * (r800 - r550))
            den = math.sqrt((2 * r800 + 1)**2 - (6 * r800 - 5 * math.sqrt(r670)) - 0.5)
            result = num / den if den!=0 else 0
            unit = "Index"

        elif formula == "wdvi":
            # [cite_start]NIR - a * Red [cite: 19]
            nir, red, a = val('NIR'), val('Red'), val('a')
            result = nir - (a * red)
            unit = "Index"

        elif formula == "pvi":
            # (NIR - a*Red - b) [cite_start]/ sqrt(a^2 + 1) [cite: 19]
            nir, red, a, b = val('NIR'), val('Red'), val('a'), val('b')
            result = (nir - a * red - b) / math.sqrt(a**2 + 1)
            unit = "Index"

        elif formula == "savi":
            nir, red, l = val('NIR'), val('Red'), val('L')
            result = ((1 + l) * (nir - red)) / (nir + red + l) if (nir+red+l)!=0 else 0
            unit = "Index"
        
        elif formula == "tsavi":
             # [cite_start]a * (R800 - a*R670 - b) / (a*R800 + R670 - a*b) [cite: 19]
             # (Note: Formula in source 19 TSAVI is complex, used simplified approximation from text)
            r800, r670, a, b = val('R800'), val('R670'), val('a'), val('b')
            num = a * (r800 - a * r670 - b)
            den = a * r800 + r670 - a * b
            result = num / den if den!=0 else 0
            unit = "Index"

        elif formula == "osavi":
             # (1 + 0.16) * (NIR - Red) [cite_start]/ (NIR + Red + 0.16) [cite: 19]
            nir, red = val('NIR'), val('Red')
            result = (1.16 * (nir - red)) / (nir + red + 0.16) if (nir+red+0.16)!=0 else 0
            unit = "Index"
        
        elif formula == "msavi":
            raw_nir = val('NIR')
            raw_red = val('Red')
            raw_L = val('L')
            print(f"DEBUG: NIR={raw_nir}, Red={raw_red}, L={raw_L}")

            # 2. Safety Check: Ensure inputs exist before doing math
            if raw_nir is None or raw_red is None:
                print("DEBUG: Error - Missing NIR or Red")
                result = 0
            else:
                try:
                    # Convert to floats
                    nir = float(raw_nir)
                    red = float(raw_red)
                    
                    # Handle L: If missing, default to 0.5. If present, convert to float.
                    if raw_L is None or raw_L == "":
                        L = 0.5
                        print("DEBUG: L was missing, defaulting to 0.5")
                    else:
                        L = float(raw_L)

                    # 3. The Math (Using the exponent 'e' from your image)
                    # Make sure 'import math' is at the top of your file!
                    term1 = (1 + L ** math.e) 
                    numerator = term1 * (nir - red)
                    denominator = nir + red + L

                    if denominator == 0:
                        result = 0
                    else:
                        result = numerator / denominator

                except Exception as e:
                    print(f"DEBUG: Calculation Error: {e}")
                    result = 0

            unit = "Index"
        
        elif formula == "msavi2":
                    # [cite_start]0.5 * (2*NIR + 1 - sqrt((2*NIR+1)^2 - 8*(NIR-Red))) [cite: 19]
             nir, red = val('NIR'), val('Red')
             term1 = 2 * nir + 1
             term2 = term1**2 - 8 * (nir - red)
             result = 0.5 * (term1 - math.sqrt(max(0, term2)))
             unit = "Index"
        
        elif formula == "sarvi":
             # Source 19 defines it with constants.
             # Assuming standard L=0.5, gamma=1
             nir, red, blue = val('R800'), val('Red'), val('Blue')
             L = 0.5
             rb = red - 1.0 * (blue - red) # Gamma approx 1
             result = ((1 + L) * (nir - rb)) / (nir + rb + L)
             unit = "Index"

        else:
            return jsonify({'error': 'Formula not implemented yet.'})

        return round(result, 4), unit
    
    except Exception as e:
        return 0, f"Error: {str(e)}"

# --- ROUTES ---

@app.route('/')
def home():
    """Renders the Single (Manual) Calculator"""
    return render_template('manual.html')

@app.route('/batch')
def batch_page():
    """Renders the Batch (Excel) Calculator"""
    return render_template('batch.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    """Endpoint for Single Manual Calculation"""
    data = request.json
    formula = data.get('formula')
    inputs = data.get('inputs')
    
    result, unit = compute_value(formula, inputs)
    
    return jsonify({'result': result, 'unit': unit})

@app.route('/upload', methods=['POST'])
def upload_file():
    """Endpoint for Batch Excel Calculation"""
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'})
    
    file = request.files['file']

    if file.filename == '':
        return jsonify({'error': 'No selected file'})

    if file:
        try:
            # 1. Read Excel
            df = pd.read_excel(file)
            
            # Convert DataFrame columns to a lowercase set for easy existence checking
            available_cols = set(col.strip().lower() for col in df.columns)
            
            calculations_performed = []

            # 2. Iterate through ALL known formulas
            for formula_name, required_cols in FORMULA_MAP.items():
                
                # Check if ALL required columns exist in the uploaded file (Case Insensitive)
                required_set = set(req.lower() for req in required_cols)
                
                if required_set.issubset(available_cols):
                    # Match found! Apply logic
                    results = df.apply(lambda row: compute_value(formula_name, row), axis=1)
                    
                    # Unzip results
                    values, units = zip(*results)
                    
                    # Add to DataFrame
                    df[f'{formula_name}_Val'] = values
                    # df[f'{formula_name}_Unit'] = units # Optional: Uncomment to see units
                    
                    calculations_performed.append(formula_name)

            if not calculations_performed:
                return jsonify({'error': 'No matching variables found. Please check column names in your Excel file against the requirements.'})

            # 3. Save result
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                df.to_excel(writer, index=False)
            output.seek(0)

            return send_file(
                output,
                as_attachment=True,
                download_name=f"Batch_Results.xlsx",
                mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
            )

        except Exception as e:
            return jsonify({'error': f"Processing failed: {str(e)}"})
    
    return jsonify({'error': 'Unknown error'})

if __name__ == '__main__':
    app.run(debug=True)
