document.addEventListener('DOMContentLoaded', () => {

    // Auth state
    let token = localStorage.getItem('token');
    let isAdmin = false;

    // Tab state
    let activeTab = 'tasks';
    let activeLaunchSubTab = 'run';

    // Tasks tab state
    let allClusterTasks = [];
    let tasksSortCol = 'state';
    let tasksSortDir = 1;
    const tasksOpenDetails = new Map();

    // Archives tab state
    let allArchives = [];
    let archivesSortCol = 'mtime';
    let archivesSortDir = -1;

    // DOM refs
    const loginContainer     = document.getElementById('login-container');
    const dashboardContainer = document.getElementById('dashboard-container');
    const loginForm          = document.getElementById('login-form');
    const loginError         = document.getElementById('login-error');
    const logoutBtn          = document.getElementById('logout-btn');

    // ── Utilities ─────────────────────────────────────────────────────────────

    function capitalize(s) {
        return s && s.length > 0
            ? s[0].toUpperCase() + s.slice(1)
            : s;
    }

    function formatTime(v) {
        if (v == null || v === '') return '';
        const secs = typeof v === 'number' ? v : parseFloat(String(v).split('/')[0]);
        if (isNaN(secs)) return '';
        const h = Math.floor(secs / 3600);
        const m = Math.floor((secs % 3600) / 60);
        const s = (secs % 60).toFixed(1);
        if (h > 0) return `${h}h${m}m${s}s`;
        if (m > 0) return `${m}m${s}s`;
        return `${s}s`;
    }

    function formatDate(mtime) {
        if (!mtime) return '';
        const d   = new Date(mtime * 1000);
        const pad = n => String(n).padStart(2, '0');
        return `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())} ${pad(d.getHours())}:${pad(d.getMinutes())}`;
    }

    function stateClass(stateLC, elevel) {
        if (stateLC === 'running')  return 'state-running';
        if (stateLC === 'queued')   return 'state-queued';
        if (stateLC === 'failed')   return 'state-finished-failed';
        if (stateLC === 'finished')
            return (elevel && elevel.startsWith('FAILED')) ? 'state-finished-failed' : 'state-finished-success';
        return '';
    }

    function stateAdjust(stateLC, elevel) {
        if (stateLC === 'finished' && elevel && elevel.startsWith('FAILED'))
            return 'failed';
        else
            return stateLC;
    }

    // ── Reusable multi-select dropdown ────────────────────────────────────────
    function buildDropdown(containerId, values, labelFn, onChange, allLabel) {
        const container = document.getElementById(containerId);
        if (!container) return null;
        container.innerHTML = '';
        const btn   = document.createElement('button');
        btn.type    = 'button';
        btn.className = 'filter-dropdown-btn';
        const panel = document.createElement('div');
        panel.className = 'dropdown-panel';
        const wrapper = document.createElement('div');
        wrapper.className = 'dropdown-wrapper';
        wrapper.appendChild(btn);
        wrapper.appendChild(panel);
        container.appendChild(wrapper);

        values.forEach(v => {
            const lbl = document.createElement('label');
            const cb  = document.createElement('input');
            cb.type = 'checkbox'; cb.className = 'dd-cb'; cb.value = v; cb.checked = true;
            cb.addEventListener('change', () => { updateBtn(); onChange(); });
            lbl.appendChild(cb);
            lbl.appendChild(document.createTextNode(' ' + (labelFn ? labelFn(v) : v)));
            panel.appendChild(lbl);
        });

        function updateBtn() {
            const all     = [...panel.querySelectorAll('.dd-cb')];
            const checked = all.filter(c => c.checked);
            btn.textContent = checked.length === all.length
                ? (allLabel || 'All') + ' \u25be'
                : `${checked.length} / ${all.length} \u25be`;
        }

        btn.addEventListener('click', e => { e.stopPropagation(); panel.classList.toggle('open'); });
        document.addEventListener('click', () => panel.classList.remove('open'));
        panel.addEventListener('click', e => e.stopPropagation());
        updateBtn();

        return {
            getSelected: () => [...panel.querySelectorAll('.dd-cb:checked')].map(c => c.value),
            resetAll: () => { panel.querySelectorAll('.dd-cb').forEach(c => c.checked = true); updateBtn(); },
        };
    }

    // ── Inline detail panel ───────────────────────────────────────────────────
    function buildDetailRow(key, content) {
        const row = document.createElement('tr');
        row.classList.add('detail-row');
        row.dataset.key = String(key);
        const td = document.createElement('td');
        td.setAttribute('colspan', '8');
        const copyBtn = document.createElement('button');
        copyBtn.textContent = 'Copy';
        copyBtn.className = 'btn-copy';
        copyBtn.addEventListener('click', () => {
            navigator.clipboard.writeText(content).then(() => {
                copyBtn.textContent = '\u2713 Copied';
                setTimeout(() => copyBtn.textContent = 'Copy', 1500);
            });
        });
        const pre = document.createElement('pre');
        pre.className = 'inline-detail';
        pre.textContent = content;
        td.appendChild(copyBtn);
        td.appendChild(pre);
        row.appendChild(td);
        return row;
    }

    function buildDetailContent(key, content) {
        const container = document.createElement('div');
        container.classList.add('detail-modal');

        // Lisible width + padding for better readability, especially on large screens
        container.style.maxWidth = '80vw';
        container.style.padding = '0px';

        const copyBtn = document.createElement('button');
        copyBtn.textContent = 'Copy';
        copyBtn.className = 'btn-copy-modal';
        copyBtn.style.marginBottom = '10px';

        copyBtn.addEventListener('click', () => {
            navigator.clipboard.writeText(content).then(() => {
                copyBtn.textContent = '\u2713 Copied';
                setTimeout(() => copyBtn.textContent = 'Copy', 1500);
            });
        });

        const pre = document.createElement('pre');
        pre.className = 'terminal';
        pre.textContent = content;

        // Scroll + word break and wrap for long content
        pre.style.maxHeight = '60vh';
        pre.style.overflow = 'auto';
        pre.style.whiteSpace = 'pre-wrap';
        pre.style.wordBreak = 'break-word';

        container.appendChild(copyBtn);
        container.appendChild(pre);

        return container;
    }

    function openModal(contentNode) {
        const modal = document.getElementById('detailModal');
        const body = document.getElementById('detailBody');

        body.innerHTML = ''; // reset
        body.appendChild(contentNode);

        modal.style.display = 'flex';
    }

    document.getElementById('detailClose').onclick = () => {
        document.getElementById('detailModal').style.display = 'none';
    };

    async function toggleDetail(openMap, row, fetchUrl, key) {
        if (openMap.has(key)) {
            openMap.delete(key);
            document.getElementById('detailModal').style.display = 'none';
            return;
        }

        try {
            const resp = await fetch(fetchUrl, { headers: { Authorization: `Bearer ${token}` } });
            const content = await resp.text();
            openMap.set(key, content);

            // 👉 ici tu réutilises ta logique existante
            openModal(buildDetailContent(key, content));

        } catch (err) {
            console.error(err);
        }
    }

    // ── Sort helpers ──────────────────────────────────────────────────────────
    function sortTasks(tasks, col, dir) {
        return [...tasks].sort((a, b) => {
            let va = a[col], vb = b[col];
            if (col === 'state') {
                const ord = { running: 0, queued: 1, waiting: 2, finished: 3, failed: 4 };
                va = ord[String(va).toLowerCase()] ?? 4;
                vb = ord[String(vb).toLowerCase()] ?? 4;
            } else if (['id', 'task_slots', 'times'].includes(col)) {
                va = Number(va) || 0; vb = Number(vb) || 0;
            } else {
                va = String(va ?? '').toLowerCase(); vb = String(vb ?? '').toLowerCase();
            }
            return va < vb ? -dir : va > vb ? dir : 0;
        });
    }

    function attachSortHeaders(tableId, getSortState, setSortState, rerender) {
        const table = document.getElementById(tableId);
        if (!table) return;
        table.querySelectorAll('th[data-col]').forEach(th => {
            th.style.cursor = 'pointer';
            th.addEventListener('click', () => {
                const col = th.dataset.col;
                const { col: cur, dir } = getSortState();
                const newDir = col === cur ? -dir : 1;
                setSortState(col, newDir);
                rerender();
                updateSortIcons(table, col, newDir);
            });
        });
    }

    function updateSortIcons(table, activeCol, dir) {
        table.querySelectorAll('th[data-col]').forEach(th => {
            const icon = th.querySelector('.sort-icon');
            if (!icon) return;
            icon.textContent = th.dataset.col === activeCol ? (dir === 1 ? ' \u25b2' : ' \u25bc') : ' \u21c5';
        });
    }

    // ── Auth ──────────────────────────────────────────────────────────────────
    async function fetchUserInfo() {
        if (!token) return;
        try {
            const resp = await fetch('/me', { headers: { Authorization: `Bearer ${token}` } });
            if (resp.ok) {
                const data = await resp.json();
                isAdmin = data.groups?.includes('admin') ?? false;
                const el = document.getElementById('current-user');
                if (el) {
                    el.textContent = data.username;
                    el.title = `Groups: ${(data.groups || []).join(', ') || 'none'}`;
                }
            }
        } catch (_) {}
    }

    function showLogin() {
        loginContainer.style.display     = 'block';
        dashboardContainer.style.display = 'none';
    }

    async function showDashboard() {
        loginContainer.style.display     = 'none';
        dashboardContainer.style.display = 'block';
        isAdmin = false;
        await fetchUserInfo();
        const launchTab = document.getElementById('tab-launch');
        if (launchTab) launchTab.style.display = isAdmin ? '' : 'none';
        switchTab('tasks');
    }

    loginForm.addEventListener('submit', async e => {
        e.preventDefault();
        const fd = new FormData();
        fd.append('username', document.getElementById('username').value);
        fd.append('password', document.getElementById('password').value);
        const resp = await fetch('/token', { method: 'POST', body: fd });
        if (resp.ok) {
            token = (await resp.json()).access_token;
            localStorage.setItem('token', token);
            await showDashboard();
        } else {
            if (loginError) loginError.textContent = 'Invalid username or password';
        }
    });

    logoutBtn.addEventListener('click', () => {
        token = null; isAdmin = false;
        localStorage.removeItem('token');
        showLogin();
    });

    if (token) { showDashboard(); } else { showLogin(); }

    // ── Tab switching ─────────────────────────────────────────────────────────
    const TAB_VIEWS = {
        tasks:    'view-tasks',
        cluster:  'view-cluster',
        launch:   'view-launch',
        archives: 'view-archives',
    };

    function switchTab(name) {
        activeTab = name;
        Object.entries(TAB_VIEWS).forEach(([tab, viewId]) => {
            const view = document.getElementById(viewId);
            const btn  = document.getElementById(`tab-${tab}`);
            if (view) view.style.display = tab === name ? '' : 'none';
            if (btn)  btn.classList.toggle('active', tab === name);
        });
        if (name === 'tasks')    fetchAllTasks();
        if (name === 'cluster')  fetchClusterSummary();
        if (name === 'archives') fetchArchives();
        if (name === 'launch')   fetchQueuesForLaunch();
    }

    async function fetchQueuesForLaunch() {
        try {
            const resp = await fetch('/queues', { headers: { Authorization: `Bearer ${token}` } });
            if (!resp.ok) return;
            const data = await resp.json();
            const dl = document.getElementById('queues-datalist');
            if (!dl) return;
            dl.innerHTML = '';
            (data.queues || []).forEach(q => {
                const opt = document.createElement('option');
                opt.value = q;
                dl.appendChild(opt);
            });
        } catch (_) { /* non-blocking */ }
    }

    ['tasks', 'cluster', 'launch', 'archives'].forEach(name => {
        const btn = document.getElementById(`tab-${name}`);
        if (btn) btn.addEventListener('click', () => switchTab(name));
    });

    // ══════════════════════════════════════════════════════════════════════════
    // TASKS TAB
    // ══════════════════════════════════════════════════════════════════════════

    let tasksNodeDd  = null;
    let tasksQueueDd = null;
    let tasksStateDd = null;

    function buildTasksFilters(tasks) {
        const nodes  = [...new Set(tasks.map(t => t.node  || ''))].filter(Boolean).sort();
        const queues = [...new Set(tasks.map(t => t.queue || ''))].filter(Boolean).sort();
        const nodesGroup  = document.getElementById('tasks-filter-nodes-group');
        const queuesGroup = document.getElementById('tasks-filter-queues-group');
        if (nodes.length > 1) {
            if (nodesGroup) nodesGroup.style.display = 'flex';
            if (!tasksNodeDd)
                tasksNodeDd = buildDropdown('tasks-filter-nodes-container', nodes, null, renderTasksTable, 'All nodes');
        }
        if (queues.length > 1) {
            if (queuesGroup) queuesGroup.style.display = 'flex';
            if (!tasksQueueDd)
                tasksQueueDd = buildDropdown('tasks-filter-queues-container', queues, null, renderTasksTable, 'All queues');
        }
        if (!tasksStateDd)
            tasksStateDd = buildDropdown('tasks-filter-states-container', ['running', 'queued', 'finished', 'failed'], null, renderTasksTable, 'All states');
    }

    function applyTasksFilters(tasks) {
        const nodes  = tasksNodeDd  ? tasksNodeDd.getSelected()  : null;
        const queues = tasksQueueDd ? tasksQueueDd.getSelected() : null;
        const states = tasksStateDd ? tasksStateDd.getSelected() : ['running', 'queued', 'finished', 'failed'];
        const search = (document.getElementById('tasks-search')?.value || '').toLowerCase();
        return tasks.filter(t => {
            if (nodes  && !nodes.includes(t.node  || ''))                      return false;
            if (queues && !queues.includes(t.queue || ''))                     return false;
            if (!states.includes((t.state || '').toLowerCase()))               return false;
            if (search && !(t.run_name || '').toLowerCase().includes(search))  return false;
            return true;
        });
    }

    async function fetchAllTasks() {
        try {
            const resp = await fetch('/cluster/tasks', { headers: { Authorization: `Bearer ${token}` } });
            if (!resp.ok) return;
            allClusterTasks = await resp.json();
            buildTasksFilters(allClusterTasks);
            renderTasksTable();
        } catch (_) {}
    }

    function renderTasksTable() {
        const tbody = document.getElementById('tasks-body');
        if (!tbody) return;
        tbody.innerHTML = '';
        const filtered = applyTasksFilters(allClusterTasks);
        const sorted   = sortTasks(filtered, tasksSortCol, tasksSortDir);
        if (!sorted.length) {
            tbody.innerHTML = '<tr><td colspan="8" style="color:#888;font-style:italic">No tasks match the current filters.</td></tr>';
            return;
        }
        sorted.forEach(task => {
            const stateLC = (task.state || '').toLowerCase();
            const cls     = stateClass(stateLC, task.elevel);
            task.state = stateAdjust(stateLC, task.elevel);
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <!-- <td title="${task.node} - ${task.node_url || ''}">${task.node || '-'}</td> -->
                <td title="${task.node} - ${task.node_url || ''}">
                    ${task.node
                        ? (task.node.length > 10 ? task.node.slice(0, 7) + '...' : task.node)
                        : '-'}
                </td>
                <td>${task.id ?? ''}</td>
                <td>${task.queue || ''}</td>
                <td class="task-slots-cell">${task.task_slots != null ? task.task_slots + ' / ' + task.queue_slots : '-'}</td>
                <td class="${cls}" title="${task.elevel || ''}">${task.state || ''}</td>
                <td>${formatTime(task.times)}</td>
                <td title="${task.run_name || ''}">
                    ${task.run_name ? (task.run_name.length > 50 ? task.run_name.slice(0, 47) + '...' : task.run_name) : '-'}
                </td>
            `;
            const actionTd = document.createElement('td');
            actionTd.className = 'action-buttons';
            const btnDefs = {
                running:  [['I', 'info'], ['L', 'log'], ['A', 'analysis'], ['K', 'kill']],
                queued:   [['I', 'info'], ['L', 'log'], ['A', 'analysis'], ['P', 'prioritize'], ['X', 'remove']],
                finished: [['I', 'info'], ['L', 'log'], ['A', 'analysis'], ['R', 'relaunch'], ['X', 'remove']],
                failed:   [['I', 'info'], ['L', 'log'], ['A', 'analysis'], ['R', 'relaunch'], ['X', 'remove']],
            };
            (btnDefs[stateLC] || [['I', 'info'], ['L', 'log'], ['A', 'analysis']]).forEach(([label, act]) => {
                const isDanger = ['kill', 'prioritize', 'remove', 'relaunch'].includes(act);
                if (isDanger && !isAdmin) return;
                const btn = document.createElement('button');
                btn.textContent = label;
                btn.title = `${capitalize(act)} '${task.run_name}'`;
                if (isDanger) btn.classList.add('btn-danger');
                if (['info', 'log', 'analysis'].includes(act)) {
                    btn.title += ` [#${task.id}] on node '${task.node || task.node_url}' and queue '${task.queue || 'default'}'`;
                    btn.addEventListener('click', () => {
                        const key    = `${task.node_url}::${task.queue}::${task.id}::${act}`;
                        const params = new URLSearchParams({ node_url: task.node_url, action: act, id: String(task.id) });
                        if (task.queue) params.append('queue', task.queue);
                        toggleDetail(tasksOpenDetails, tr, `/cluster/proxy/queue?${params}`, key);
                    });
                } else {
                    if (act === 'relaunch') {
                        btn.title += ` on queue '${task.queue || 'default'}'`;
                    } else {
                        btn.title += ` [#${task.id}] on node '${task.node || task.node_url}' and queue '${task.queue || 'default'}'`;
                    }
                    btn.addEventListener('click', async () => {
                        if (!confirm(`${capitalize(act)} '${task.run_name}' [#${task.id}] on queue '${task.queue || 'default'}' on node '${task.node || task.node_url}'?`)) return;
                        btn.disabled = true; btn.textContent = '\u2026';
                        let ok = false;
                        try {
                            let resp;
                            if (act === 'relaunch' || act === 'prioritize') {
                                const p = new URLSearchParams({ node_url: task.node_url });
                                if (task.queue) p.append('queue', task.queue);
                                if (act === 'prioritize') p.append('prioritize', true);
                                resp = await fetch(`/cluster/proxy/relaunch/${task.id}?${p}`, {
                                    method: 'POST',
                                    headers: { Authorization: `Bearer ${token}` },
                                });
                            } else {
                                const p = new URLSearchParams({ node_url: task.node_url, action: act, id: String(task.id) });
                                if (task.queue) p.append('queue', task.queue);
                                resp = await fetch(`/cluster/proxy/queue?${p}`, {
                                    headers: { Authorization: `Bearer ${token}` },
                                });
                            }
                            ok = resp.ok;
                        } catch (_) {}
                        // btn.textContent = ok ? '\u2713 Done' : '\u2717 Failed';
                        btn.textContent = ok ? '\u2713 Done' : '\u2717 Failed';
                        btn.classList.toggle('btn-success', ok);
                        btn.classList.toggle('btn-error',   !ok);
                        setTimeout(() => {
                            btn.textContent = label; btn.disabled = false;
                            btn.classList.remove('btn-success', 'btn-error');
                            fetchAllTasks();
                        }, 1500);
                    });
                }
                actionTd.appendChild(btn);
            });
            tr.appendChild(actionTd);
            tbody.appendChild(tr);
        });
    }

    attachSortHeaders(
        'tasks-table',
        () => ({ col: tasksSortCol, dir: tasksSortDir }),
        (col, dir) => { tasksSortCol = col; tasksSortDir = dir; },
        renderTasksTable
    );

    document.getElementById('tasks-filter-reset')?.addEventListener('click', () => {
        tasksNodeDd?.resetAll(); tasksQueueDd?.resetAll(); tasksStateDd?.resetAll();
        const s = document.getElementById('tasks-search'); if (s) s.value = '';
        renderTasksTable();
    });
    document.getElementById('tasks-search')?.addEventListener('input', renderTasksTable);

    // ══════════════════════════════════════════════════════════════════════════
    // CLUSTER TAB
    // ══════════════════════════════════════════════════════════════════════════

    async function fetchClusterSummary() {
        const el = document.getElementById('cluster-summary-content');
        if (!el) return;
        try {
            const resp = await fetch('/cluster/summary', { headers: { Authorization: `Bearer ${token}` } });
            if (!resp.ok) { el.innerHTML = `<p class="error">Error ${resp.status}</p>`; return; }
            renderClusterSummary(await resp.json(), el);
        } catch (e) { el.innerHTML = `<p class="error">Could not reach /cluster/summary: ${e}</p>`; }
    }

    function renderClusterSummary(data, el) {
        const nodes  = data.nodes  || {};
        const totals = data.totals || {};
        if (!Object.keys(nodes).length) {
            el.innerHTML = '<div class="cluster-no-peers">No peers configured. Add entries to <code>config/peers.json</code> to enable cluster view.</div>';
            return;
        }
        const table = document.createElement('table');
        table.className = 'cluster-table';
        table.innerHTML = `<thead><tr>
            <th>Node</th><th>Queue</th>
            <th title="Total slots">Config</th>
            <th title="Running">Running</th>
            <th title="Queued">Queued</th>
            <th title="Available">Available</th>
            <th title="Usage">Usage</th>
        </tr></thead>`;
        const tbody = document.createElement('tbody');
        for (const [name, node] of Object.entries(nodes)) {
            const urlStr = node.url
                ? `<span style="font-weight:400;font-size:0.82em;color:#888;margin-left:6px">${node.url}</span>` : '';
            if (node.status !== 'online') {
                const tr = document.createElement('tr');
                tr.className = 'cluster-node-offline';
                tr.innerHTML = `<td colspan="7">${name}<span class="node-status offline">offline</span><br>${urlStr}</td>`;
                tbody.appendChild(tr); continue;
            }
            const queues      = node.queues || {};
            const queueNames  = Object.keys(queues);
            queueNames.forEach((qname, i) => {
                const q  = queues[qname];
                const tr = document.createElement('tr');
                if (i === 0) {
                    const nodeTd = document.createElement('td');
                    nodeTd.rowSpan = queueNames.length;
                    nodeTd.innerHTML = `<strong>${name}</strong><span class="node-status online">online</span><br>${urlStr}`;
                    tr.appendChild(nodeTd);
                }
                const used   = (q.running || 0) + (q.queued || 0);
                const pct    = q.configured > 0 ? Math.round(used / q.configured * 100) : 0;
                const barCls = pct >= 100 ? 'full' : pct >= 70 ? 'warn' : '';
                tr.innerHTML += `<td>${qname}</td><td>${q.configured}</td>
                    <td class="state-running">${q.running}</td>
                    <td class="state-queued">${q.queued}</td>
                    <td>${q.available}</td>
                    <td><div class="slot-bar-wrap">
                        <span class="slot-bar"><span class="slot-bar-fill ${barCls}" style="width:${Math.min(pct, 100)}%"></span></span>
                        <span class="slot-num">${pct}%</span></div></td>`;
                tbody.appendChild(tr);
            });
        }
        for (const [qname, t] of Object.entries(totals)) {
            const used   = (t.running || 0) + (t.queued || 0);
            const pct    = t.configured > 0 ? Math.round(used / t.configured * 100) : 0;
            const barCls = pct >= 100 ? 'full' : pct >= 70 ? 'warn' : '';
            const tr = document.createElement('tr');
            tr.style.cssText = 'font-weight:700;background:#f5f6f7;border-top:2px solid #aaa';
            tr.innerHTML = `<td>TOTAL</td><td>${qname}</td><td>${t.configured}</td>
                <td class="state-running">${t.running}</td>
                <td class="state-queued">${t.queued}</td>
                <td>${t.available}</td>
                <td><div class="slot-bar-wrap">
                    <span class="slot-bar"><span class="slot-bar-fill ${barCls}" style="width:${Math.min(pct, 100)}%"></span></span>
                    <span class="slot-num">${pct}%</span></div></td>`;
            tbody.appendChild(tr);
        }
        table.appendChild(tbody);
        el.innerHTML = ''; el.appendChild(table);
    }

    // ══════════════════════════════════════════════════════════════════════════
    // LAUNCH TAB
    // ══════════════════════════════════════════════════════════════════════════

    document.getElementById('launch-tab-run')?.addEventListener('click', () => {
        activeLaunchSubTab = 'run';
        document.getElementById('launch-view-run').style.display     = '';
        document.getElementById('launch-view-advanced').style.display = 'none';
        document.getElementById('launch-view-command-docker').style.display = 'none';
        document.getElementById('launch-tab-run').classList.add('active');
        document.getElementById('launch-tab-advanced').classList.remove('active');
        document.getElementById('launch-tab-command-docker').classList.remove('active');
        document.getElementById('launch-result').style.display = 'none';
    });
    document.getElementById('launch-tab-advanced')?.addEventListener('click', () => {
        activeLaunchSubTab = 'advanced';
        document.getElementById('launch-view-run').style.display     = 'none';
        document.getElementById('launch-view-advanced').style.display = '';
        document.getElementById('launch-view-command-docker').style.display = 'none';
        document.getElementById('launch-tab-run').classList.remove('active');
        document.getElementById('launch-tab-advanced').classList.add('active');
        document.getElementById('launch-tab-command-docker').classList.remove('active');
        document.getElementById('launch-result').style.display = 'none';
    });
    document.getElementById('launch-tab-command-docker')?.addEventListener('click', () => {
        activeLaunchSubTab = 'command-docker';
        document.getElementById('launch-view-run').style.display     = 'none';
        document.getElementById('launch-view-advanced').style.display = 'none';
        document.getElementById('launch-view-command-docker').style.display = '';
        document.getElementById('launch-tab-run').classList.remove('active');
        document.getElementById('launch-tab-advanced').classList.remove('active');
        document.getElementById('launch-tab-command-docker').classList.add('active');
        document.getElementById('launch-result').style.display = 'none';
    });

    async function submitAnalysis(payload) {
        const resp = await fetch('/analysis', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${token}` },
            body: payload,
        });
        const text = await resp.text();
        showLaunchResult(text, resp.ok);
        if (resp.ok) setTimeout(fetchAllTasks, 500);
    }

    function showLaunchResult(text, success) {
        const box     = document.getElementById('launch-result');
        const content = document.getElementById('launch-result-content');
        if (!box || !content) return;
        box.style.display = '';
        const m = text.match(/^(STARK\.\w+\.(ID-\w+-NAME-(.+)))$/);
        if (m && success) {
            content.innerHTML = `<div class="launch-result-ok">
                <div class="launch-result-row"><span class="lr-label">Analysis ID</span><code>${m[1]}</code></div>
                <div class="launch-result-row"><span class="lr-label">Name</span><span>${m[3] || '-'}</span></div>
            </div>`;
        } else {
            content.innerHTML = `<div class="${success ? 'launch-result-ok' : 'launch-result-err'}">${text}</div>`;
        }
    }

    document.getElementById('analysis-form-run')?.addEventListener('submit', async e => {
        e.preventDefault();
        const runName = document.getElementById('run-name-input')?.value.trim();
        if (!runName) return;
        await submitAnalysis(JSON.stringify({ run: runName }));
    });

    document.getElementById('analysis-form-advanced')?.addEventListener('submit', async e => {
        e.preventDefault();
        const payload = document.getElementById('json-input')?.value.trim();
        if (!payload) return;
        try { JSON.parse(payload); } catch (_) { showLaunchResult('Invalid JSON payload', false); return; }
        await submitAnalysis(payload);
    });

    document.getElementById('analysis-form-command-docker')?.addEventListener('submit', async e => {
        e.preventDefault();
        const analysis_name = document.getElementById('command-docker-analysis-name')?.value.trim();
        if (!analysis_name) {
            showLaunchResult('Analysis name is required', false);
            return;
        }
        const image = document.getElementById('command-docker-image')?.value.trim();
        if (!analysis_name || !image) {
            showLaunchResult('Docker image is required', false);
            return;
        }
        const command = document.getElementById('command-docker-command-docker')?.value.trim();
        if (!command) {
            showLaunchResult('Docker command is required', false);
            return;
        }
        const docker_extra_params = document.getElementById('command-docker-extra-params')?.value.trim();
        const use_stark_container_mount = document.getElementById('command-docker-use-stark-container-mount')?.checked;
        const queue = document.getElementById('command-docker-queue')?.value.trim();
        const threads = Number.parseInt(document.getElementById('command-docker-threads')?.value.trim());
        const memory = document.getElementById('command-docker-memory')?.value.trim();
        const prioritize = document.getElementById('command-docker-prioritize')?.checked;
        await submitAnalysis(JSON.stringify({ "analysis_name": analysis_name, "image": image, "command_docker": command, "docker_extra_params": docker_extra_params, "use_stark_container_mount": use_stark_container_mount, "queue": queue, "threads": threads, "memory": memory, "prioritize": prioritize }));
    });

    // ══════════════════════════════════════════════════════════════════════════
    // ARCHIVES TAB
    // ══════════════════════════════════════════════════════════════════════════

    let archivesNodeDd   = null;
    let archivesQueueDd  = null;
    let archivesStatusDd = null;

    function buildArchivesFilters(archives) {
        const nodes    = [...new Set(archives.map(a => a.node   || ''))].filter(Boolean).sort();
        const queues   = [...new Set(archives.map(a => a.queue  || ''))].filter(Boolean).sort();
        const statuses = [...new Set(archives.map(a => a.status || 'unknown'))].sort();
        const nodesGroup  = document.getElementById('archives-filter-nodes-group');
        const queuesGroup = document.getElementById('archives-filter-queues-group');
        if (nodes.length > 1) {
            if (nodesGroup) nodesGroup.style.display = 'flex';
            if (!archivesNodeDd)
                archivesNodeDd = buildDropdown('archives-filter-nodes-container', nodes, null, renderArchivesTable, 'All nodes');
        }
        if (queues.length > 1) {
            if (queuesGroup) queuesGroup.style.display = 'flex';
            if (!archivesQueueDd)
                archivesQueueDd = buildDropdown('archives-filter-queues-container', queues, null, renderArchivesTable, 'All queues');
        }
        if (!archivesStatusDd)
            archivesStatusDd = buildDropdown('archives-filter-status-container', statuses, null, renderArchivesTable, 'All statuses');
    }

    function applyArchivesFilters(archives) {
        const nodes    = archivesNodeDd   ? archivesNodeDd.getSelected()   : null;
        const queues   = archivesQueueDd  ? archivesQueueDd.getSelected()  : null;
        const statuses = archivesStatusDd ? archivesStatusDd.getSelected() : null;
        const search   = (document.getElementById('archives-search')?.value || '').toLowerCase();
        return archives.filter(a => {
            if (nodes    && !nodes.includes(a.node   || ''))           return false;
            if (queues   && !queues.includes(a.queue  || ''))           return false;
            if (statuses && !statuses.includes(a.status || 'unknown')) return false;
            if (search   && !(a.run_name || '').toLowerCase().includes(search)) return false;
            return true;
        });
    }

    function sortByCol(arr, col, dir) {
        return [...arr].sort((a, b) => {
            let va = a[col], vb = b[col];
            if (['mtime', 'threads'].includes(col)) { va = Number(va) || 0; vb = Number(vb) || 0; }
            else { va = String(va ?? '').toLowerCase(); vb = String(vb ?? '').toLowerCase(); }
            return va < vb ? -dir : va > vb ? dir : 0;
        });
    }

    function renderArchivesTable() {
        const el = document.getElementById('archives-content');
        if (!el) return;
        const filtered = applyArchivesFilters(allArchives);
        const sorted   = sortByCol(filtered, archivesSortCol, archivesSortDir);
        if (!sorted.length) {
            el.innerHTML = '<p style="color:#888;font-style:italic">No archives found.</p>'; return;
        }
        // const statusCls = { done: 'state-finished-success', failed: 'state-finished-failed', unknown: 'state-unknown' };
        const statusCls = { finished: 'state-finished-success', failed: 'state-finished-failed', unknown: 'state-unknown' };
        const table = document.createElement('table');
        table.className = 'cluster-tasks-table';
        table.innerHTML = `<thead><tr>
            <th data-col="queue">Queue <span class="sort-icon"></span></th>
            <th data-col="threads">Slots <span class="sort-icon"></span></th>
            <th data-col="status">Status <span class="sort-icon"></span></th>
            <th data-col="mtime">Date <span class="sort-icon"></span></th>
            <th data-col="run_name">Analysis Name <span class="sort-icon"></span></th>
        </tr></thead>`;
        updateSortIcons(table, archivesSortCol, archivesSortDir);
        table.querySelectorAll('th[data-col]').forEach(th => {
            th.style.cursor = 'pointer';
            th.addEventListener('click', () => {
                const col = th.dataset.col;
                archivesSortDir = col === archivesSortCol ? -archivesSortDir : -1;
                archivesSortCol = col;
                renderArchivesTable();
            });
        });
        const tbody = document.createElement('tbody');
        sorted.forEach(a => {
            const tr = document.createElement('tr');
            const statusLabel = a.status || 'unknown';
            let date_title = a.launch_date ? `Launch:\t${formatDate(a.launch_date) ?? '-'}` : 'Date unknown';
            date_title += a.end_date ? `\nEnd:\t\t${formatDate(a.end_date) ?? '-'}` : '';
            date_title += a.exec_time ? `\nTime:\t${formatTime(a.exec_time) ?? '-'}` : '';
            tr.innerHTML = `
                <td>${a.queue || '-'}</td>
                <td>${a.threads ?? '-'}</td>
                <td class="${statusCls[statusLabel] || 'state-unknown'}">${statusLabel}</td>
                <td title="${date_title}">${formatDate(a.mtime) ?? '-'}</td>
                <td title="${a.analysis_id_name || '-'}" >${a.run_name || '-'}</td>
            `;
            tbody.appendChild(tr);
        });
        table.appendChild(tbody);
        el.innerHTML = ''; el.appendChild(table);
    }

    // Attach event listeners for archives tab
    document.getElementById('tasks-refresh-btn')?.addEventListener('click', fetchAllTasks);

    // Cluster summary is less volatile, so no auto-refresh, only manual
    document.getElementById('cluster-refresh-btn')?.addEventListener('click', fetchClusterSummary);

    // Archives can be heavy to load, so no auto-refresh, only manual
    document.getElementById('archives-refresh-btn')?.addEventListener('click', fetchArchives);
    document.getElementById('archives-filter-reset')?.addEventListener('click', () => {
        archivesNodeDd?.resetAll(); archivesQueueDd?.resetAll(); archivesStatusDd?.resetAll();
        const s = document.getElementById('archives-search'); if (s) s.value = '';
        renderArchivesTable();
    });
    document.getElementById('archives-search')?.addEventListener('input', renderArchivesTable);

    async function fetchArchives() {
        const el = document.getElementById('archives-content');
        if (!el) return;
        el.innerHTML = '<p class="cluster-loading">Loading\u2026</p>';
        try {
            const resp = await fetch('/cluster/archives', { headers: { Authorization: `Bearer ${token}` } });
            if (!resp.ok) { el.innerHTML = `<p class="error">Error ${resp.status}</p>`; return; }
            allArchives = await resp.json();
            // uniquify json response by node+queue+analysis_id_name to avoid duplicates from multiple nodes
            const seen = new Set();
            allArchives = allArchives.filter(a => {
                const key = `${a.analysis_id_name}`;
                if (seen.has(key)) return false;
                seen.add(key);
                return true;
            });
            archivesNodeDd   = null;
            archivesQueueDd  = null;
            archivesStatusDd = null;
            buildArchivesFilters(allArchives);
            renderArchivesTable();
        } catch (e) { el.innerHTML = `<p class="error">Could not reach /cluster/archives: ${e}</p>`; }
    }

    // ── Auto-refresh ──────────────────────────────────────────────────────────
    const REFRESH_MS = typeof REFRESH_INTERVAL_MS !== 'undefined' ? REFRESH_INTERVAL_MS : 10000;
    setInterval(() => {
        if (activeTab === 'tasks')   fetchAllTasks();
        if (activeTab === 'cluster') fetchClusterSummary();
    }, REFRESH_MS);

});
