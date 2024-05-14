async function sendPostRequest(url, body, headers = {'Content-Type': 'application/x-www-form-urlencoded'}) {
    const response = await fetch(url, {
        method: 'POST',
        headers: headers,
        body: body
    });
    const data = await response.json();
    if (response.ok) {
        return data;
    } else {
        throw new Error(data.detail || 'Unknown error occurred');
    }
}

async function executeQuery() {
    const databaseName = document.getElementById('database-name').value;
    const query = document.getElementById('any-query').value;
    const body = `database_name=${encodeURIComponent(databaseName)}&query=${encodeURIComponent(query)}`;
    const url = 'http://localhost:8000/execute-query/';

    try {
        const result = await sendPostRequest(url, body);
        document.getElementById('result-query').textContent = JSON.stringify(result);
    } catch (error) {
        document.getElementById('result-query').textContent = 'Error: ' + error.message;
    }
}