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
    const openDetails = new Map();
    let allTasks = []; // full unfiltered task list

    // ── Filters ────────────────────────────────────────────────────────────────
    function getActiveStates() {
        return [...document.querySelectorAll('.filter-state:checked')].map(el => el.value);
    }
    function getActiveQueues() {
        const checks = document.querySelectorAll('.filter-queue:checked');
        if (!checks.length) return null; // no queue checkboxes yet → no filter
        return [...checks].map(el => el.value);
    }
    function isFailedOnly() {
        return document.getElementById('filter-failed-only')?.checked ?? false;
    }

    function applyFilters() {
        const states = getActiveStates();
        const queues = getActiveQueues();
        const failed = isFailedOnly();
        return allTasks.filter(task => {
            const stateLC = task.state.toLowerCase();
            // 'Failed only' = show exclusively finished+failed tasks
            // if (failed) return stateLC === 'finished' && task.elevel && task.elevel.startsWith('FAILED');
            if (!states.includes(stateLC)) return false;
            if (queues && !queues.includes(task.queue || '')) return false;
            // 'Failed only' applies only to finished tasks; running/queued always pass
            if (failed && stateLC === 'finished' && !(task.elevel && task.elevel.startsWith('FAILED'))) return false;
            return true;
        });
    }

    function buildStateDropdown() {
        const container = document.getElementById('filter-states-container');
        if (!container || document.getElementById('state-dropdown')) return;
        const items = [
            { value: 'running',  label: 'Running',  cls: 'filter-state', checked: true },
            { value: 'queued',   label: 'Queued',   cls: 'filter-state', checked: true },
            { value: 'finished', label: 'Finished', cls: 'filter-state', checked: true },
            null,
            { value: 'failed-only', label: 'Failed only (finished)', id: 'filter-failed-only', cls: 'filter-failed-only-cb', checked: false },
        ];
        const btn = document.createElement('button');
        btn.id = 'state-dropdown-btn';
        btn.type = 'button';
        btn.className = 'filter-dropdown-btn';
        const panel = document.createElement('div');
        panel.id = 'state-dropdown-panel';
        panel.className = 'dropdown-panel';
        items.forEach(item => {
            if (item === null) {
                const hr = document.createElement('hr');
                hr.className = 'dropdown-divider';
                panel.appendChild(hr);
                return;
            }
            const lbl = document.createElement('label');
            const cb  = document.createElement('input');
            cb.type = 'checkbox'; cb.className = item.cls; cb.value = item.value; cb.checked = item.checked;
            if (item.id) cb.id = item.id;
            cb.addEventListener('change', () => { updateStateDropdownBtn(); renderFilteredTable(); });
            lbl.appendChild(cb);
            lbl.appendChild(document.createTextNode(' ' + item.label));
            panel.appendChild(lbl);
        });
        const wrapper = document.createElement('div');
        wrapper.id = 'state-dropdown';
        wrapper.className = 'dropdown-wrapper';
        wrapper.appendChild(btn);
        wrapper.appendChild(panel);
        container.appendChild(wrapper);
        btn.addEventListener('click', (e) => { e.stopPropagation(); panel.classList.toggle('open'); });
        document.addEventListener('click', () => panel.classList.remove('open'));
        panel.addEventListener('click', e => e.stopPropagation());
        updateStateDropdownBtn();
    }

    function updateStateDropdownBtn() {
        const btn = document.getElementById('state-dropdown-btn');
        if (!btn) return;
        const states  = [...document.querySelectorAll('.filter-state')];
        const checked = states.filter(cb => cb.checked);
        const failed  = document.getElementById('filter-failed-only')?.checked;
        let label = checked.length === states.length ? 'All states' : `${checked.length} / ${states.length} states`;
        if (failed) label += ' · Failed only';
        btn.textContent = label + ' ▾';
    }

    function rebuildQueueCheckboxes(tasks) {
        const group = document.getElementById('filter-queues-group');
        if (!group) return;
        const seen = new Set(tasks.map(t => t.queue || ''));
        if (seen.size <= 1) { group.style.display = 'none'; return; }
        group.style.display = 'flex';

        // Build dropdown once, then only add new queues
        let dropdown = document.getElementById('queue-dropdown');
        if (!dropdown) {
            // Button that toggles the panel
            const btn = document.createElement('button');
            btn.id = 'queue-dropdown-btn';
            btn.type = 'button';
            // Panel
            const panel = document.createElement('div');
            panel.id = 'queue-dropdown-panel';
            panel.className = 'dropdown-panel';
            // Wrapper
            dropdown = document.createElement('div');
            dropdown.id = 'queue-dropdown';
            dropdown.className = 'dropdown-wrapper';
            dropdown.appendChild(btn);
            dropdown.appendChild(panel);
            const container = document.getElementById('filter-queues-checkboxes');
            if (container) container.appendChild(dropdown);

            btn.addEventListener('click', (e) => {
                e.stopPropagation();
                panel.classList.toggle('open');
            });
            document.addEventListener('click', () => panel.classList.remove('open'));
            panel.addEventListener('click', e => e.stopPropagation());
        }

        const panel = document.getElementById('queue-dropdown-panel');
        const existing = new Set([...panel.querySelectorAll('.filter-queue')].map(el => el.value));
        seen.forEach(q => {
            if (existing.has(q)) return;
            const lbl = document.createElement('label');
            const cb  = document.createElement('input');
            cb.type = 'checkbox'; cb.className = 'filter-queue'; cb.value = q; cb.checked = true;
            cb.addEventListener('change', () => { updateDropdownBtn(); renderFilteredTable(); });
            lbl.appendChild(cb);
            lbl.appendChild(document.createTextNode(' ' + (q || '(default)')));
            panel.appendChild(lbl);
        });
        updateDropdownBtn();
    }

    function updateDropdownBtn() {
        const btn = document.getElementById('queue-dropdown-btn');
        if (!btn) return;
        const all   = [...document.querySelectorAll('.filter-queue')];
        const checked = all.filter(cb => cb.checked);
        btn.textContent = checked.length === all.length
            ? 'All queues ▾'
            : `${checked.length} / ${all.length} queues ▾`;
    }

    function renderFilteredTable() {
        const queueBody = document.getElementById('queue-body');
        if (!queueBody) return;
        const filtered = applyFilters();
        queueBody.innerHTML = '';
        if (filtered.length > 0) {
            document.getElementById('queue-table-container').style.display = 'block';
            filtered.forEach(task => {
                const row = document.createElement('tr');
                const stateLC = task.state.toLowerCase();
                let stateClass = '';
                if (stateLC === 'running') stateClass = 'state-running';
                else if (stateLC === 'queued') stateClass = 'state-queued';
                else if (stateLC === 'finished') stateClass = task.elevel && task.elevel.startsWith('FAILED') ? 'state-finished-failed' : 'state-finished-success';
                row.innerHTML = `
                    <td>${task.id}</td>
                    <td class="${stateClass}">${task.state}</td>
                    <td>${task.queue || ''}</td>
                    <td class="task-slots-cell" title="${task.task_slots != null ? 'Uses ' + task.task_slots + ' / ' + task.queue_slots + ' slots' : 'Slot info unavailable'}">${task.task_slots != null ? task.task_slots + ' / ' + task.queue_slots : '\u2014'}</td>
                    <td>${task.elevel}</td>
                    <td>${formatTime(task.times)}</td>
                    <td>${task.run_name}</td>
                `;
                const actionTd = document.createElement('td');
                actionTd.className = 'action-buttons';
                const stateButtons = {
                    'running':  [['I', 'info'], ['L', 'log'], ['A', 'analysis'], ['K', 'kill']],
                    'queued':   [['I', 'info'], ['L', 'log'], ['A', 'analysis'], ['P', 'prioritize'], ['X', 'remove']],
                    'finished': [['I', 'info'], ['L', 'log'], ['A', 'analysis'], ['R', 'relaunch'], ['X', 'remove']],
                };
                const buttons = stateButtons[task.state.toLowerCase()] || [['I', 'info'], ['L', 'log'], ['A', 'analysis']];
                buttons.forEach(([label, act]) => {
                    if (['kill', 'prioritize', 'remove', 'relaunch'].includes(act) && !isAdmin) return;
                    const btn = document.createElement('button');
                    btn.textContent = label;
                    btn.title = `${label} for task ${task.id} on queue ${task.queue || 'default'}`;
                    if (['kill', 'prioritize', 'remove', 'relaunch'].includes(act)) {
                        btn.classList.add('btn-danger');
                    }
                    const actionLabels = { kill: 'Kill', prioritize: 'Prioritize', remove: 'Remove', relaunch: 'Relaunch' };
                    if (['info', 'analysis', 'log'].includes(act)) {
                        btn.addEventListener('click', () => toggleInlineDetail(row, act, task.id, task.queue));
                    } else if (['kill', 'prioritize', 'remove', 'relaunch'].includes(act)) {
                        btn.addEventListener('click', async () => {
                            const actionName = actionLabels[act] || act;
                            if (!confirm(`${actionName} task #${task.id} on queue "${task.queue || 'default'}" ?`)) return;
                            btn.disabled = true;
                            const originalText = btn.textContent;
                            btn.textContent = '…';
                            let ok = false;
                            try {
                                let resp;
                                if (act === 'relaunch') {
                                    const qParam = task.queue ? `?queue=${encodeURIComponent(task.queue)}` : '';
                                    resp = await fetch(`/relaunch/${task.id}${qParam}`, {
                                        method: 'POST',
                                        headers: { 'Authorization': `Bearer ${token}` }
                                    });
                                } else {
                                    let url = `/queue?action=${act}&id=${task.id}`;
                                    if (task.queue) url += `&queue=${encodeURIComponent(task.queue)}`;
                                    resp = await fetch(url, { headers: { 'Authorization': `Bearer ${token}` } });
                                }
                                ok = resp.ok;
                            } catch (_) { ok = false; }
                            btn.textContent = ok ? '\u2713 Done' : '\u2717 Failed';
                            btn.classList.toggle('btn-success', ok);
                            btn.classList.toggle('btn-error', !ok);
                            setTimeout(() => {
                                btn.textContent = originalText;
                                btn.disabled = false;
                                btn.classList.remove('btn-success', 'btn-error');
                                getQueue('list');
                            }, 1500);
                        });
                    } else {
                        btn.addEventListener('click', () => getQueue(act, task.id, task.queue));
                    }
                    actionTd.appendChild(btn);
                });
                row.appendChild(actionTd);
                queueBody.appendChild(row);
            });
        } else {
            queueBody.innerHTML = '<tr><td colspan="7">No tasks match the current filters.</td></tr>';
        }
        restoreOpenDetails();
    }

    // ── Filter event wiring ────────────────────────────────────────────────────
    document.querySelectorAll('.filter-state').forEach(cb => cb.addEventListener('change', renderFilteredTable));
    document.getElementById('filter-failed-only')?.addEventListener('change', renderFilteredTable);
    document.getElementById('filter-reset')?.addEventListener('click', () => {
        document.querySelectorAll('.filter-state').forEach(cb => cb.checked = true);
        document.querySelectorAll('.filter-queue').forEach(cb => cb.checked = true);
        const fo = document.getElementById('filter-failed-only');
        if (fo) fo.checked = false;
        updateStateDropdownBtn();
        updateDropdownBtn();
        renderFilteredTable();
    });

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

    async function fetchUserInfo() {
        if (!token) return;
        const response = await fetch('/me', {
            headers: { 'Authorization': `Bearer ${token}` }
        });
        if (response.ok) {
            const data = await response.json();
            isAdmin = data.groups && data.groups.includes('admin');
            const el = document.getElementById('current-user');
            if (el) {
                el.textContent = data.username;
                el.title = `Groups: ${(data.groups || []).join(', ') || 'none'}`;
            }
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
        isAdmin = false;
        await fetchUserInfo();
        const launchSection = document.getElementById('analysis-form-container');
        if (launchSection) launchSection.style.display = isAdmin ? 'block' : 'none';
        buildStateDropdown();
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

    async function getQueue(action = 'list', id = '', queue = '') {

        // For list action, we fetch and build the table
        if (action === 'list') {
            const response = await fetch('/queue?action=list', {
                headers: { 'Authorization': `Bearer ${token}` }
            });
            const queueOutput = document.getElementById('queue-output');

            if (response.ok) {
                const tasks = await response.json();
                // const stateOrder = { 'running': 0, 'queued': 1, 'finished': 2 };
                // tasks.sort((a, b) => {
                //     const oa = stateOrder[a.state.toLowerCase()] ?? 3;
                //     const ob = stateOrder[b.state.toLowerCase()] ?? 3;
                //     return oa - ob;
                // });
                const stateOrder = { 'running': 0, 'queued': 1, 'finished': 2 };
                tasks.sort((a, b) => {
                    const oa = stateOrder[a.state.toLowerCase()] ?? 3;
                    const ob = stateOrder[b.state.toLowerCase()] ?? 3;
                    if (oa !== ob) return oa - ob;
                    return Number(a.id) - Number(b.id);
                });
                allTasks = tasks;
                rebuildQueueCheckboxes(tasks);
                renderFilteredTable();
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
            if (id) url += `&id=${id}`;
            if (queue) url += `&queue=${encodeURIComponent(queue)}`;
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
        detailTd.setAttribute('colspan', '7');

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
        openDetails.forEach((detail, key) => {
            const sep = key.indexOf('::');
            const keyQueue = key.substring(0, sep);
            const keyId = key.substring(sep + 2);
            let found = false;
            for (const row of taskRows) {
                const idCell = row.querySelector('td:first-child');
                const queueCell = row.querySelector('td:nth-child(3)');
                if (idCell && idCell.textContent.trim() === keyId &&
                    queueCell && queueCell.textContent.trim() === keyQueue) {
                    const next = row.nextElementSibling;
                    if (!next || !next.classList.contains('detail-row')) {
                        row.insertAdjacentElement('afterend', buildDetailRow(key, detail.content));
                    }
                    found = true;
                    break;
                }
            }
            if (!found) toDelete.push(key);
        });
        toDelete.forEach(k => openDetails.delete(k));
    }

    async function toggleInlineDetail(row, action, id, queue = '') {
        const key = `${queue}::${id}`;
        const existingDetail = row.nextElementSibling;
        if (existingDetail && existingDetail.classList.contains('detail-row')) {
            existingDetail.remove();
            openDetails.delete(key);
            return;
        }

        const queueParam = queue ? `&queue=${encodeURIComponent(queue)}` : '';
        const url = `/queue?action=${action}&id=${id}${queueParam}`;
        const response = await fetch(url, {
            headers: { 'Authorization': `Bearer ${token}` }
        });
        const content = await response.text();

        openDetails.set(key, { action, queue, content });
        row.insertAdjacentElement('afterend', buildDetailRow(key, content));
    }

    // Auto-refresh the queue every 10 seconds
    setInterval(() => getQueue('list'), typeof REFRESH_INTERVAL_MS !== 'undefined' ? REFRESH_INTERVAL_MS : 10000);
});