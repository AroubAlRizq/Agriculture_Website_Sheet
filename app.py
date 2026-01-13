from flask import Flask, render_template, request, jsonify, send_file
import math
import pandas as pd
import io

app = Flask(__name__)

# --- CORE CALCULATION ENGINE ---
def compute_value(formula, inputs):
    """
    Computes the result for a single row (dictionary-like inputs).
    Returns a tuple: (result_value, unit_string)
    """
    result = 0
    unit = ""

    # Helper: Safe Float Conversion from the row
    # We strip whitespace from keys just in case the Excel header has spaces "NIR "
    def val(key):
        try:
            # Case insensitive lookup attempt
            for k in inputs.keys():
                if k.strip().lower() == key.lower():
                    return float(inputs[k])
            return 0.0
        except:
            return 0.0

    try:
        # ==========================================
        # 1. EVAPOTRANSPIRATION (ET)
        # ==========================================
        if formula == "hargreaves":
            t_mean, t_max, t_min, ra = val('T_mean'), val('T_max'), val('T_min'), val('Ra')
            result = 0.0023 * (t_mean + 17.8) * math.sqrt(max(0, t_max - t_min)) * ra
            unit = "mm/day"
        
        elif formula == "blaney_criddle":
            ta, kc = val('Ta'), val('Kc')
            result = (0.0173 * ta - 0.314) * kc * ta
            unit = "mm/day"

        elif formula == "fao56":
            rn, g, t, u2 = val('Rn'), val('G'), val('T_mean'), val('u2')
            es_ea, delta, gamma = val('es_ea'), val('Delta'), val('Gamma')
            num = (0.408 * delta * (rn - g)) + (gamma * (900 / (t + 273)) * u2 * es_ea)
            den = delta + (gamma * (1 + 0.34 * u2))
            result = num / den if den != 0 else 0
            unit = "mm/day"

        elif formula == "stephens_stewart":
            ta, rl = val('Ta'), val('Rl')
            result = (0.0082 * ta - 0.19) * (rl / 1500.0) * 25.4
            unit = "mm/day"

        elif formula == "grassi":
            ta, rl = val('Ta'), val('Rl')
            result = 0.537 * 0.000675 * rl * (0.62 + 0.00559 * ta) * 25.4
            unit = "mm/day"

        elif formula == "linarce":
            ta, td, z, lat = val('Ta'), val('Td'), val('z'), val('lat')
            tm = ta + 0.006 * z
            num = ((700 * tm) / (100 - lat)) + (15 * (ta - td))
            den = 80 - ta
            result = num / den if den != 0 else 0
            unit = "mm/day"

        # ==========================================
        # 2. GROWING DEGREE DAYS (GDD)
        # ==========================================
        elif formula == "gdd_arnold":
            tm_big, tm_small, tb = val('TM'), val('Tm'), val('Tb')
            result = ((tm_big + tm_small) / 2) - tb
            if result < 0: result = 0
            unit = "deg C-days"

        elif formula == "gdd_villa_nova":
            TM, Tm, Tb = val('TM'), val('Tm'), val('Tb')
            if Tb >= TM:
                result = 0
            elif Tb < Tm:
                result = ((Tm - Tb) + (TM - Tm)) / 2.0 
            elif Tb >= Tm: 
                result = ((TM - Tb)**2) / (2 * (TM - Tm)) if (TM-Tm) != 0 else 0
            else: result = 0
            unit = "deg C-days"

        elif formula == "gdd_ometto":
            TM, Tm, Tb, TB = val('TM'), val('Tm'), val('Tb'), val('TB')
            if TM > TB and TB > Tm and Tm > Tb: # Case 4
                den = 2 * (TM - Tm)
                result = (2 * (TM - Tm) * (Tm - Tb) + (TM - Tm)**2 - (TM - TB)**2) / den if den != 0 else 0
            elif TB > TM and TM > Tm and Tm > Tb: # Case 1
                result = ((TM - Tm) / 2.0) + (Tm - Tb)
            elif TM > TB and TB > Tb and Tb > Tm: # Case 5
                den = 2 * (TM - Tm)
                result = 0.5 * (((TM - Tb)**2 - (TM - TB)**2) / (TM - Tm)) if (TM-Tm) != 0 else 0
            elif TB > TM and TM > Tb and Tb > Tm: # Case 2
                 result = ((TM - Tb)**2) / (2 * (TM - Tm)) if (TM-Tm) != 0 else 0
            else:
                result = 0
            unit = "deg C-days"

        elif formula == "gdd_snyder":
            TM, Tm, Tb, TB = val('TM'), val('Tm'), val('Tb'), val('TB')
            M = (TM + Tm) / 2.0
            W = (TM - Tm) / 2.0
            
            def get_angle(thresh):
                if W == 0: return 0
                val_c = (thresh - M) / W
                val_c = max(-1, min(1, val_c))
                return math.asin(val_c)

            theta = get_angle(Tb)
            phi = get_angle(TB)
            pi = math.pi

            if Tb <= Tm: gdd1 = M - Tb
            elif Tb >= TM: gdd1 = 0
            else: gdd1 = ((1/pi) * ((M - Tb) * (pi/2 - theta) + (W * math.cos(theta))))

            if TB <= Tm: gdd2 = M - TB
            elif TB >= TM: gdd2 = 0
            else: gdd2 = ((1/pi) * ((M - TB) * (pi/2 - phi) + (W * math.cos(phi))))

            result = gdd1 - gdd2
            unit = "deg C-days"

        # ==========================================
        # 3. ACCUMULATION CHILLING UNITS
        # ==========================================
        elif formula == "chill_utah":
            t = val('T_current')
            if t <= 1.4: result = 0
            elif 1.5 <= t <= 2.4: result = 0.5
            elif 2.5 <= t <= 9.1: result = 1.0
            elif 9.2 <= t <= 12.4: result = 0.5
            elif 12.5 <= t <= 15.9: result = 0
            elif 16.0 <= t <= 18.0: result = -0.5
            elif t > 18.0: result = -1.0
            unit = "Chill Units"

        elif formula == "chill_nc":
            t = val('T_current')
            if t <= 1.5: result = 0
            elif 1.6 <= t <= 7.1: result = 0.5
            elif 7.2 <= t <= 12.9: result = 1.0
            elif 13.0 <= t <= 16.4: result = 0.5
            elif 16.5 <= t <= 19.0: result = 0
            elif 19.1 <= t <= 20.6: result = -0.5
            elif 20.7 <= t <= 22.0: result = -1.0
            elif 22.1 <= t <= 23.2: result = -1.5
            else: result = -2.0
            unit = "Chill Units"

        # ==========================================
        # 4. VEGETATION INDICES
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
            nir, red_edge = val('R790'), val('R720')
            result = (nir - red_edge) / (nir + red_edge) if (nir+red_edge)!=0 else 0
            unit = "Index"

        elif formula == "ccci":
            ndre, n_min, n_max = val('NDRE'), val('NDRE_min'), val('NDRE_max')
            result = (ndre - n_min) / (n_max - n_min) if (n_max-n_min)!=0 else 0
            unit = "Index"

        elif formula == "rvi":
            nir, red = val('NIR'), val('Red')
            result = nir / red if red!=0 else 0
            unit = "Ratio"

        elif formula == "evi":
            nir, red, blue = val('NIR'), val('Red'), val('Blue')
            den = nir + (6 * red) - (7.5 * blue) + 1
            result = 2.5 * ((nir - red) / den) if den!=0 else 0
            unit = "Index"

        elif formula == "evi2":
            nir, red = val('NIR'), val('Red')
            den = nir + (2.4 * red) + 1
            result = 2.5 * ((nir - red) / den) if den!=0 else 0
            unit = "Index"
        
        elif formula == "varigreen":
            green, red, blue = val('Green'), val('Red'), val('Blue')
            den = green + red - blue
            result = (green - red) / den if den!=0 else 0
            unit = "Index"
        
        elif formula == "vari700":
            r700, red, blue = val('R700'), val('Red'), val('Blue')
            num = r700 - (1.7 * red) + (0.7 * blue)
            den = r700 + (2.3 * red) - (1.3 * blue)
            result = num / den if den!=0 else 0
            unit = "Index"

        elif formula == "tvi":
            r750, r550, r670 = val('R750'), val('R550'), val('R670')
            result = 0.5 * (120 * (r750 - r550) - 200 * (r670 - r550))
            unit = "Index"

        elif formula == "mtvi1":
            r800, r550, r670 = val('R800'), val('R550'), val('R670')
            result = 1.2 * (1.2 * (r800 - r550) - 2.5 * (r670 - r550))
            unit = "Index"
            
        elif formula == "mtvi2":
            r800, r550, r670 = val('R800'), val('R550'), val('R670')
            num = 1.5 * (1.2 * (r800 - r550) - 2.5 * (r670 - r550))
            den = math.sqrt(max(0, (2 * r800 + 1)**2 - (6 * r800 - 5 * math.sqrt(max(0, r670))) - 0.5))
            result = num / den if den!=0 else 0
            unit = "Index"

        elif formula == "mtci":
            r753, r708, r681 = val('R753'), val('R708'), val('R681')
            den = r708 - r681
            result = (r753 - r708) / den if den!=0 else 0
            unit = "Index"
        
        elif formula == "car":
            r700, r670, r550 = val('R700'), val('R670'), val('R550')
            a = (r700 - r550) / 150.0 
            b = r550 - (a * 550)    
            term1 = abs(a * 670 + b + r670)
            term2 = math.sqrt(a**2 + 1)
            result = term1/term2 if term2!=0 else 0
            unit = "Index"

        elif formula == "cari":
            r700, r670, r550 = val('R700'), val('R670'), val('R550')
            a = (r700 - r550) / 150.0 
            b = r550 - (a * 550)    
            term1 = abs(a * 670 + b + r670)
            term2 = math.sqrt(a**2 + 1)
            term3 = term1/term2 if term2!=0 else 0
            term4 = r700/r670 if r670!=0 else 0
            result = term3*term4
            unit = "Index"

        elif formula == "tcari":
            r700, r670, r550 = val('R700'), val('R670'), val('R550')
            term1 = r700 - r670
            term2 = 0.2 * (r700 - r550) * (r700 / r670) if r670!=0 else 0
            result = 3 * (term1 - term2)
            unit = "Index"

        elif formula == "mcari":
            r700, r670, r550 = val('R700'), val('R670'), val('R550')
            result = ((r700 - r670) - 0.2 * (r700 - r550)) * (r700 / r670) if r670!=0 else 0
            unit = "Index"

        elif formula == "mcari1":
            r800, r670, r550 = val('R800'), val('R670'), val('R550')
            result = 1.2 * (2.5 * (r800 - r670) - 1.3 * (r800 - r550))
            unit = "Index"

        elif formula == "mcari2":
            r800, r670, r550 = val('R800'), val('R670'), val('R550')            
            num = 1.5 * (2.5 * (r800 - r670) - 1.3 * (r800 - r550))
            den = math.sqrt(max(0, (2 * r800 + 1)**2 - (6 * r800 - 5 * math.sqrt(max(0, r670))) - 0.5))
            result = num / den if den!=0 else 0
            unit = "Index"

        elif formula == "wdvi":
            nir, red, a = val('NIR'), val('Red'), val('a')
            result = nir - (a * red)
            unit = "Index"

        elif formula == "pvi":
            nir, red, a, b = val('NIR'), val('Red'), val('a'), val('b')
            result = (nir - a * red - b) / math.sqrt(a**2 + 1)
            unit = "Index"

        elif formula == "savi":
            nir, red, l = val('NIR'), val('Red'), val('L')
            if l == 0: l = 0.5 # Default fallback
            result = ((1 + l) * (nir - red)) / (nir + red + l) if (nir+red+l)!=0 else 0
            unit = "Index"
        
        elif formula == "tsavi":
            r800, r670, a, b = val('R800'), val('R670'), val('a'), val('b')
            num = a * (r800 - a * r670 - b)
            den = a * r800 + r670 - a * b
            result = num / den if den!=0 else 0
            unit = "Index"

        elif formula == "osavi":
            nir, red = val('NIR'), val('Red')
            result = (1.16 * (nir - red)) / (nir + red + 0.16) if (nir+red+0.16)!=0 else 0
            unit = "Index"
        
        elif formula == "msavi":
            nir, red, L = val('NIR'), val('Red'), val('L')
            if L == 0: L = 0.5
            term1 = (1 + L ** math.e) 
            num = term1 * (nir - red)
            den = nir + red + L
            result = num / den if den != 0 else 0
            unit = "Index"
        
        elif formula == "msavi2":
             nir, red = val('NIR'), val('Red')
             term1 = 2 * nir + 1
             term2 = term1**2 - 8 * (nir - red)
             result = 0.5 * (term1 - math.sqrt(max(0, term2)))
             unit = "Index"
        
        elif formula == "sarvi":
             nir, red, blue = val('R800'), val('Red'), val('Blue')
             L = 0.5
             rb = red - 1.0 * (blue - red)
             result = ((1 + L) * (nir - rb)) / (nir + rb + L)
             unit = "Index"

        else:
            return None, "Error"

        return round(result, 4), unit
    
    except Exception as e:
        return 0, f"Error: {str(e)}"

# --- ROUTES ---

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'})
    
    file = request.files['file']
    formula = request.form.get('formula')

    if file.filename == '':
        return jsonify({'error': 'No selected file'})

    if file and formula:
        try:
            # 1. Read Excel
            df = pd.read_excel(file)
            
            # 2. Apply Computation Row by Row
            # We use a lambda to call compute_value for each row
            results = df.apply(lambda row: compute_value(formula, row), axis=1)
            
            # 3. Create new columns
            # unzip the results into two lists
            values, units = zip(*results)
            
            df[f'{formula}_Result'] = values
            df[f'{formula}_Unit'] = units

            # 4. Save back to Excel in memory
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                df.to_excel(writer, index=False)
            output.seek(0)

            # 5. Return the file
            return send_file(
                output,
                as_attachment=True,
                download_name=f"calculated_{formula}.xlsx",
                mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
            )

        except Exception as e:
            return jsonify({'error': f"Processing failed: {str(e)}"})
    
    return jsonify({'error': 'Unknown error'})

if __name__ == '__main__':
    app.run(debug=True)
