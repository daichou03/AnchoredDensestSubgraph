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