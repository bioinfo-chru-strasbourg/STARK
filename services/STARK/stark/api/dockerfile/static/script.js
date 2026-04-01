document.addEventListener('DOMContentLoaded', () => {
    const loginContainer = document.getElementById('login-container');
    const dashboardContainer = document.getElementById('dashboard-container');
    const loginForm = document.getElementById('login-form');
    const loginError = document.getElementById('login-error');
    const analysisForm = document.getElementById('analysis-form');
    const analysisResponse = document.getElementById('analysis-response');
    const queueOutput = document.getElementById('queue-output');
    const logoutBtn = document.getElementById('logout-btn');

    let token = localStorage.getItem('token');
    let isAdmin = false;
    const openDetails = new Map(); // taskId -> { action, content }

    async function fetchUserInfo() {
        if (!token) return;
        const response = await fetch('/me', {
            headers: { 'Authorization': `Bearer ${token}` }
        });
        if (response.ok) {
            const data = await response.json();
            isAdmin = data.groups && data.groups.includes('admin');
        }
    }

    if (token) {
        showDashboard();
    } else {
        showLogin();
    }

    function showLogin() {
        loginContainer.style.display = 'block';
        dashboardContainer.style.display = 'none';
    }

    async function showDashboard() {
        loginContainer.style.display = 'none';
        dashboardContainer.style.display = 'block';
        isAdmin = false; // reset before fetching
        await fetchUserInfo();
        const launchSection = document.getElementById('analysis-form-container');
        if (launchSection) launchSection.style.display = isAdmin ? 'block' : 'none';
        getQueue();
    }

    loginForm.addEventListener('submit', async (e) => {
        e.preventDefault();
        const username = document.getElementById('username').value;
        const password = document.getElementById('password').value;

        const formData = new FormData();
        formData.append('username', username);
        formData.append('password', password);

        const response = await fetch('/token', {
            method: 'POST',
            body: formData
        });

        if (response.ok) {
            const data = await response.json();
            token = data.access_token;
            localStorage.setItem('token', token);
            showDashboard();
        } else {
            loginError.textContent = 'Invalid username or password';
        }
    });

    logoutBtn.addEventListener('click', () => {
        token = null;
        isAdmin = false;
        localStorage.removeItem('token');
        showLogin();
    });

    analysisForm.addEventListener('submit', async (e) => {
        e.preventDefault();
        const runName = document.getElementById('run-name-input').value.trim();
        const jsonInput = document.getElementById('json-input').value.trim();

        let payload;
        if (jsonInput) {
            payload = jsonInput;
        } else if (runName) {
            payload = JSON.stringify({ run: runName });
        } else {
            analysisResponse.textContent = 'Please enter a run name or a JSON payload.';
            return;
        }
        const response = await fetch(`/analysis`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${token}`
            },
            body: payload
        });
        const data = await response.text();
        analysisResponse.textContent = data;
        getQueue();
    });

    async function getQueue(action = 'list', id = '') {

        function formatTime(timesStr) {
            if (!timesStr) return '';
            if (timesStr === 'N/A') return 'N/A';
            const firstTime = parseFloat(timesStr.split('/')[0]);
            if (isNaN(firstTime)) return 'N/A';
            const h = Math.floor(firstTime / 3600);
            const m = Math.floor((firstTime % 3600) / 60);
            const s = (firstTime % 60).toFixed(1);
            if (h > 0) return `${h}h${m}m${s}s`;
            if (m > 0) return `${m}m${s}s`;
            return `${s}s`;
        }

        // For list action, we fetch and build the table
        if (action === 'list') {
            const response = await fetch('/queue?action=list', {
                headers: { 'Authorization': `Bearer ${token}` }
            });
            const queueBody = document.getElementById('queue-body');
            const queueOutput = document.getElementById('queue-output');
            
            if (response.ok) {
                const tasks = await response.json();
                queueBody.innerHTML = ''; // Clear table
                if (tasks.length > 0) {
                    document.getElementById('queue-table-container').style.display = 'block';
                    tasks.forEach(task => {
                        const row = document.createElement('tr');
                        row.innerHTML = `
                            <td>${task.id}</td>
                            <td>${task.state}</td>
                            <td>${task.elevel}</td>
                            <td>${formatTime(task.times)}</td>
                            <td>${task.run_name}</td>
                        `;
                        const actionTd = document.createElement('td');
                        actionTd.className = 'action-buttons';
                        const stateButtons = {
                            'running':  [['Info', 'info'], ['Log', 'log'], ['Analysis', 'analysis'], ['Kill', 'kill']],
                            'queued':   [['Info', 'info'], ['Log', 'log'], ['Analysis', 'analysis'], ['Prioritize', 'prioritize'], ['Remove', 'remove']],
                            'finished': [['Info', 'info'], ['Log', 'log'], ['Analysis', 'analysis'], ['Relaunch', 'relaunch']],
                        };
                        const buttons = stateButtons[task.state.toLowerCase()] || [['Info', 'info'], ['Log', 'log'], ['Analysis', 'analysis']];
                        buttons.forEach(([label, act]) => {
                            if (['kill', 'prioritize', 'remove', 'relaunch'].includes(act) && !isAdmin) return;
                            const btn = document.createElement('button');
                            btn.textContent = label;
                            if (['kill', 'prioritize', 'remove', 'relaunch'].includes(act)) {
                                btn.classList.add('btn-danger');
                            }
                            if (['info', 'analysis', 'log'].includes(act)) {
                                btn.addEventListener('click', () => toggleInlineDetail(row, act, task.id));
                            } else if (act === 'relaunch') {
                                btn.addEventListener('click', async () => {
                                    btn.disabled = true;
                                    btn.textContent = 'Launching...';
                                    const resp = await fetch(`/relaunch/${task.id}`, {
                                        method: 'POST',
                                        headers: { 'Authorization': `Bearer ${token}` }
                                    });
                                    const text = await resp.text();
                                    btn.textContent = resp.ok ? '\u2713 Done' : '\u2717 Failed';
                                    if (resp.ok) setTimeout(() => getQueue('list'), 1500);
                                });
                            } else {
                                btn.addEventListener('click', () => getQueue(act, task.id));
                            }
                            actionTd.appendChild(btn);
                        });
                        row.appendChild(actionTd);
                        queueBody.appendChild(row);
                    });
                } else {
                    queueBody.innerHTML = '<tr><td colspan="6">No tasks in the queue.</td></tr>';
                }
                restoreOpenDetails();
            } else {
                // Handle errors
                const errorText = await response.text();
                document.getElementById('queue-table-container').style.display = 'none';
                queueOutput.style.display = 'block';
                queueOutput.textContent = `Error fetching queue: ${errorText}`;
            }
        } else {
            // For other actions (kill, remove, prioritize)
            let url = `/queue?action=${action}`;
            if (id) {
                url += `&id=${id}`;
            }
            const response = await fetch(url, {
                headers: { 'Authorization': `Bearer ${token}` }
            });
            await response.text();

            // Only refresh the list for state-changing actions
            if (['kill', 'remove', 'prioritize'].includes(action)) {
                setTimeout(() => getQueue('list'), 1000);
            }
        }
    }

    function buildDetailRow(taskId, content) {
        const detailRow = document.createElement('tr');
        detailRow.classList.add('detail-row');
        detailRow.dataset.taskId = String(taskId);
        const detailTd = document.createElement('td');
        detailTd.setAttribute('colspan', '6');

        const copyBtn = document.createElement('button');
        copyBtn.textContent = 'Copy';
        copyBtn.className = 'btn-copy';
        copyBtn.addEventListener('click', () => {
            navigator.clipboard.writeText(content).then(() => {
                copyBtn.textContent = '✓ Copied';
                setTimeout(() => copyBtn.textContent = 'Copy', 1500);
            });
        });

        const pre = document.createElement('pre');
        pre.className = 'inline-detail';
        pre.textContent = content;

        detailTd.appendChild(copyBtn);
        detailTd.appendChild(pre);
        detailRow.appendChild(detailTd);
        return detailRow;
    }

    function restoreOpenDetails() {
        const taskRows = document.querySelectorAll('#queue-body tr:not(.detail-row)');
        const toDelete = [];
        openDetails.forEach((detail, taskId) => {
            let found = false;
            for (const row of taskRows) {
                const idCell = row.querySelector('td:first-child');
                if (idCell && idCell.textContent.trim() === String(taskId)) {
                    const next = row.nextElementSibling;
                    if (!next || !next.classList.contains('detail-row')) {
                        row.insertAdjacentElement('afterend', buildDetailRow(taskId, detail.content));
                    }
                    found = true;
                    break;
                }
            }
            if (!found) toDelete.push(taskId);
        });
        toDelete.forEach(id => openDetails.delete(id));
    }

    async function toggleInlineDetail(row, action, id) {
        const taskId = String(id);
        const existingDetail = row.nextElementSibling;
        if (existingDetail && existingDetail.classList.contains('detail-row')) {
            existingDetail.remove();
            openDetails.delete(taskId);
            return;
        }

        const url = `/queue?action=${action}&id=${id}`;
        const response = await fetch(url, {
            headers: { 'Authorization': `Bearer ${token}` }
        });
        const content = await response.text();

        openDetails.set(taskId, { action, content });
        row.insertAdjacentElement('afterend', buildDetailRow(taskId, content));
    }

    // Auto-refresh the queue every 10 seconds
    setInterval(() => getQueue('list'), typeof REFRESH_INTERVAL_MS !== 'undefined' ? REFRESH_INTERVAL_MS : 10000);
});