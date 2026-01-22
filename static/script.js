/* =========================================
   CONFIGURATION: Inputs required for Manual Mode
   (Matches Python FORMULA_MAP)
   ========================================= */
const FORMULA_INPUTS = {
    // ET
    'hargreaves': ['T_mean', 'T_max', 'T_min', 'Ra'],
    'blaney_criddle': ['Ta', 'Kc'],
    'fao56': ['Rn', 'G', 'T_mean', 'u2', 'es_ea', 'Delta', 'Gamma'],
    'stephens_stewart': ['Ta', 'Rl'],
    'grassi': ['Ta', 'Rl'],
    'linarce': ['Ta', 'Td', 'z', 'lat'],

    // GDD
    'gdd_arnold': ['TM', 'Tm', 'Tb'],
    'gdd_villa_nova': ['TM', 'Tm', 'Tb'],
    'gdd_ometto': ['TM', 'Tm', 'Tb', 'TB'],
    'gdd_snyder': ['TM', 'Tm', 'Tb', 'TB'],

    // Chill
    'chill_utah': ['T_current'],
    'chill_nc': ['T_current'],

    // Vegetation Indices (All require specific bands)
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
};

document.addEventListener('DOMContentLoaded', function() {
    console.log("Agricultural Calculator Script Loaded.");

    // =========================================
    // LOGIC 1: MANUAL CALCULATOR PAGE
    // =========================================
    const formulaSelect = document.getElementById('formula-select');
    const inputsContainer = document.getElementById('inputs-container');
    const manualCalcBtn = document.getElementById('calc-btn'); // Note: ID is same on both pages
    
    // Check if we are on the Manual Page
    if (formulaSelect && inputsContainer) {
        
        // 1. Handle Dropdown Change (Generate Input Fields)
        formulaSelect.addEventListener('change', function() {
            const formula = this.value;
            inputsContainer.innerHTML = ''; // Clear previous
            document.getElementById('result-area').classList.add('hidden');

            const requiredInputs = FORMULA_INPUTS[formula];

            if (requiredInputs) {
                requiredInputs.forEach(key => {
                    const html = `
                        <div class="input-bubble">
                            <label>${key}</label>
                            <input type="number" step="any" name="${key}" placeholder="0.00">
                        </div>
                    `;
                    inputsContainer.insertAdjacentHTML('beforeend', html);
                });
            }
        });

        // 2. Handle Calculate Button Click
        manualCalcBtn.addEventListener('click', async function() {
            const formula = formulaSelect.value;
            if (!formula) {
                alert("Please select a formula first.");
                return;
            }

            // Gather Inputs
            const inputs = {};
            const inputFields = inputsContainer.querySelectorAll('input');
            inputFields.forEach(field => {
                inputs[field.name] = field.value;
            });

            // Send to Server
            try {
                const response = await fetch('/calculate', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ formula: formula, inputs: inputs })
                });

                const data = await response.json();

                // Show Result
                const resArea = document.getElementById('result-area');
                const resVal = document.getElementById('result-display');
                const resUnit = document.getElementById('unit-display');

                resArea.classList.remove('hidden');
                resVal.innerText = data.result;
                resUnit.innerText = data.unit;

            } catch (error) {
                console.error("Error:", error);
                alert("Calculation failed. Check console.");
            }
        });
    }

    // =========================================
    // LOGIC 2: BATCH PROCESSING PAGE
    // =========================================
    const fileInput = document.getElementById('excel-file');
    const loadingMsg = document.getElementById('loading-msg');

    // Check if we are on the Batch Page (Excel input exists)
    if (fileInput) {
        
        // Reuse manualCalcBtn variable because the ID is 'calc-btn' on both pages
        manualCalcBtn.addEventListener('click', async function() {
            const file = fileInput.files[0];

            if (!file) {
                alert("Please upload an Excel file first.");
                return;
            }

            // UI Feedback
            const originalText = manualCalcBtn.innerText;
            manualCalcBtn.disabled = true;
            manualCalcBtn.style.opacity = "0.7";
            manualCalcBtn.innerText = "Processing...";
            if(loadingMsg) loadingMsg.style.display = "block";

            const formData = new FormData();
            formData.append('file', file);

            try {
                const response = await fetch('/upload', {
                    method: 'POST',
                    body: formData
                });

                if (response.ok) {
                    const blob = await response.blob();
                    const url = window.URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.style.display = 'none';
                    a.href = url;
                    a.download = `Batch_Results_${new Date().toISOString().slice(0,10)}.xlsx`; 
                    document.body.appendChild(a);
                    a.click();
                    window.URL.revokeObjectURL(url);
                    
                    alert("Success! Your processed file has been downloaded.");
                } else {
                    const errorData = await response.json();
                    alert("Error: " + (errorData.error || "Unknown server error"));
                }

            } catch (error) {
                console.error('Error:', error);
                alert("Failed to communicate with the server.");
            } finally {
                manualCalcBtn.disabled = false;
                manualCalcBtn.style.opacity = "1";
                manualCalcBtn.innerText = originalText;
                if(loadingMsg) loadingMsg.style.display = "none";
            }
        });
    }
});
