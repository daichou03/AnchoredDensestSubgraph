export async function sendPostRequest(url, body, headers = {'Content-Type': 'application/x-www-form-urlencoded'}) {
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