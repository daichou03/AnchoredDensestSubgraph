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

async function getNodeCount() {
    const databaseName = document.getElementById('database-name').value;
    const body = `database_name=${encodeURIComponent(databaseName)}`;
    const url = 'http://localhost:8000/get-node-count/';

    try {
        const data = await sendPostRequest(url, body);
        document.getElementById('result-node-count').textContent = 'Node count in "' + databaseName + '": ' + data.node_count;
    } catch (error) {
        document.getElementById('result-node-count').textContent = 'Error: ' + error.message;
    }
}