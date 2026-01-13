document.addEventListener('DOMContentLoaded', function() {
    console.log("Agricultural Calculator Script Loaded.");

    const fileInput = document.getElementById('excel-file');
    const calcBtn = document.getElementById('calc-btn');
    const loadingMsg = document.getElementById('loading-msg');

    calcBtn.addEventListener('click', async function() {
        const file = fileInput.files[0];

        if (!file) {
            alert("Please upload an Excel file first.");
            return;
        }

        // UI Feedback
        calcBtn.disabled = true;
        calcBtn.style.opacity = "0.7";
        calcBtn.innerText = "Processing...";
        loadingMsg.style.display = "block";

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
                a.download = `Smart_Results_${new Date().toISOString().slice(0,10)}.xlsx`; 
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
                
                alert("Success! All possible indices have been calculated and downloaded.");
            } else {
                const errorData = await response.json();
                alert("Error: " + (errorData.error || "Unknown server error"));
            }

        } catch (error) {
            console.error('Error:', error);
            alert("Failed to communicate with the server.");
        } finally {
            calcBtn.disabled = false;
            calcBtn.style.opacity = "1";
            calcBtn.innerText = "Auto-Process Data";
            loadingMsg.style.display = "none";
        }
    });
});
