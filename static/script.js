document.addEventListener('DOMContentLoaded', function() {
    console.log("Agricultural Calculator Script Loaded.");

    const select = document.getElementById('formula-select');
    const fileInput = document.getElementById('excel-file');
    const calcBtn = document.getElementById('calc-btn');
    const loadingMsg = document.getElementById('loading-msg');

    calcBtn.addEventListener('click', async function() {
        // Validation
        const formula = select.value;
        const file = fileInput.files[0];

        if (!formula) {
            alert("Please select a Mathematical Model.");
            return;
        }
        if (!file) {
            alert("Please upload an Excel file.");
            return;
        }

        // UI Feedback
        calcBtn.disabled = true;
        calcBtn.style.opacity = "0.7";
        calcBtn.innerText = "Computing...";
        loadingMsg.style.display = "block";

        // Prepare Data
        const formData = new FormData();
        formData.append('formula', formula);
        formData.append('file', file);

        try {
            const response = await fetch('/upload', {
                method: 'POST',
                body: formData
            });

            if (response.ok) {
                // Handle File Download (Blob)
                const blob = await response.blob();
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.style.display = 'none';
                a.href = url;
                // Name the file based on the formula
                a.download = `Results_${formula}.xlsx`; 
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
                
                // Reset UI
                alert("Computation Complete! File downloaded.");
            } else {
                const errorData = await response.json();
                alert("Error: " + (errorData.error || "Unknown server error"));
            }

        } catch (error) {
            console.error('Error:', error);
            alert("Failed to communicate with the server.");
        } finally {
            // Restore UI
            calcBtn.disabled = false;
            calcBtn.style.opacity = "1";
            calcBtn.innerText = "Process & Download Results";
            loadingMsg.style.display = "none";
        }
    });
});
