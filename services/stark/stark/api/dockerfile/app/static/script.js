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
    let _currentOpenDetail = null; // { map, key } — tracks which openMap entry is visible in the modal

    // Archives tab state
    let allArchives = [];
    const archivesOpenDetails = new Map();
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
        if (_currentOpenDetail) {
            _currentOpenDetail.map.delete(_currentOpenDetail.key);
            _currentOpenDetail = null;
        }
    };

    async function toggleDetail(openMap, row, fetchUrl, key) {
        if (openMap.has(key)) {
            openMap.delete(key);
            _currentOpenDetail = null;
            document.getElementById('detailModal').style.display = 'none';
            return;
        }

        // Close any previously open detail from a different button
        if (_currentOpenDetail) {
            _currentOpenDetail.map.delete(_currentOpenDetail.key);
            _currentOpenDetail = null;
        }

        try {
            const resp = await fetch(fetchUrl, { headers: { Authorization: `Bearer ${token}` } });
            const content = await resp.text();
            openMap.set(key, content);
            _currentOpenDetail = { map: openMap, key };
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
        const launchEnabled = typeof LAUNCH_ENABLED !== 'undefined' ? LAUNCH_ENABLED : true;
        const launchModes   = typeof LAUNCH_MODES   !== 'undefined' ? LAUNCH_MODES   : ['run', 'analysis', 'docker', 'advanced'];
        const launchTab = document.getElementById('tab-launch');
        if (launchTab) launchTab.style.display = (isAdmin && launchEnabled) ? '' : 'none';
        // Apply sub-tab visibility according to LAUNCH_MODES
        const modeButtonMap = { run: 'launch-tab-run', analysis: 'launch-tab-analysis', docker: 'launch-tab-docker', advanced: 'launch-tab-advanced' };
        Object.entries(modeButtonMap).forEach(([mode, btnId]) => {
            const btn = document.getElementById(btnId);
            if (btn) btn.style.display = launchModes.includes(mode) ? '' : 'none';
        });
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
        stats:    'view-stats',
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
        if (name === 'stats')    fetchStats();
    }

    async function fetchQueuesForLaunch() {
        // Ensure the active sub-tab is one of the allowed modes; if not, switch to the first allowed
        const launchModes = typeof LAUNCH_MODES !== 'undefined' ? LAUNCH_MODES : ['run', 'analysis', 'docker', 'advanced'];
        if (!launchModes.includes(activeLaunchSubTab)) {
            const firstMode = launchModes[0];
            document.getElementById(`launch-tab-${firstMode || 'run'}`)?.click();
        }
        try {
            const resp = await fetch('/queues', { headers: { Authorization: `Bearer ${token}` } });
            if (!resp.ok) return;
            const data = await resp.json();
            // Legacy datalist for docker form
            const dl = document.getElementById('queues-datalist');
            if (dl) {
                dl.innerHTML = '';
                (data.queues || []).forEach(q => {
                    const opt = document.createElement('option');
                    opt.value = q;
                    dl.appendChild(opt);
                });
            }
            // STARK Analysis queue uses an <input list="queues-datalist">,
            // so its suggestions must come from the datalist populated above.
        } catch (_) { /* non-blocking */ }

        // Load modules for STARK Analysis form
        try {
            const mresp = await fetch('/modules', { headers: { Authorization: `Bearer ${token}` } });
            if (!mresp.ok) return;
            const mdata = await mresp.json();
            const minput = document.getElementById('analysis-module');
            const mdl = document.getElementById('modules-datalist');
            if (minput && mdl) {
                const previousValue = minput.value;
                const alreadyLoaded = !!window._modulesMeta;
                mdl.innerHTML = '';
                window._modulesMeta = {};
                (mdata.modules || []).forEach(m => {
                    if (m.enable === false || m.available === false) return; // skip disabled or unavailable modules
                    const opt = document.createElement('option');
                    opt.value = m.name;
                    opt.label = m.description ? `${m.description}` : m.name;
                    mdl.appendChild(opt);
                    window._modulesMeta[m.name] = { defaults: m.defaults || {} };
                });
                // Restore previous value if still valid, otherwise leave empty
                const modules = mdata.modules || [];
                const stillValid = previousValue && modules.some(m => m.name === previousValue);
                if (stillValid) minput.value = previousValue;
                // Apply defaults only on first load (not on tab switch)
                if (!alreadyLoaded) applyAnalysisModuleDefaults();
                if (!alreadyLoaded) minput.addEventListener('input', applyAnalysisModuleDefaults);
            }
        } catch (_) { /* non-blocking */ }
    }

    ['tasks', 'cluster', 'launch', 'archives', 'stats'].forEach(name => {
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
            const analysis_name_label = task.run_name ? (task.run_name.length > 50 ? task.run_name.slice(0, 47) + '...' : task.run_name) : '-';
            const analysis_title = task.run_name ? `Analysis Name:\t${task.run_name}` : '';
            const node_title = task.node ? `Node:\t${task.node}\nURL:\t\t${task.node_url || '-'}` : 'Node unknown';
            const node_label = task.node ? (task.node.length > 12 ? task.node.slice(0, 12) + '...' : task.node) : '-';
            tr.innerHTML = `
                <!-- <td title="${task.node} - ${task.node_url || ''}">${task.node || '-'}</td> -->
                <td title="${node_title}">${node_label}</td>
                <td>${task.id ?? ''}</td>
                <td>${task.queue || ''}</td>
                <td class="task-slots-cell">${task.task_slots != null ? task.task_slots + ' / ' + task.queue_slots : '-'}</td>
                <td class="${cls}" title="${task.elevel || ''}">${task.state || ''}</td>
                <td>${formatTime(task.times)}</td>
                <td title="${analysis_title}" >${analysis_name_label}</td>
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
                        btn.disabled = true; btn.textContent = '...';
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
            el.innerHTML = '<div class="cluster-no-nodes">No nodes configured. Add entries to <code>config/nodes.json</code> to enable cluster view.</div>';
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

    // Analysis name: character filter + counter (shared logic for both forms)
    const ANALYSIS_NAME_MAX = 80;
    const ANALYSIS_NAME_RE  = /[^A-Za-z0-9._-]/g;
    function setupAnalysisNameInput(inputId, hintId) {
        const input = document.getElementById(inputId);
        const hint  = document.getElementById(hintId);
        if (!input || !hint) return;
        function update() {
            const filtered = input.value.replace(ANALYSIS_NAME_RE, '_');
            if (filtered !== input.value) {
                const pos = input.selectionStart;
                input.value = filtered;
                input.setSelectionRange(pos, pos);
            }
            const len = input.value.length;
            hint.textContent = `${len} / ${ANALYSIS_NAME_MAX}`;
            hint.classList.toggle('input-hint-warn', len >= ANALYSIS_NAME_MAX * 0.9);
        }
        input.addEventListener('input', update);
        input.addEventListener('paste', () => setTimeout(update, 0));
        update();
    }
    setupAnalysisNameInput('analysis-name-input',  'analysis-name-input-hint');
    setupAnalysisNameInput('docker-analysis-name', 'docker-analysis-name-hint');

    document.getElementById('launch-tab-run')?.addEventListener('click', () => {
        activeLaunchSubTab = 'run';
        document.getElementById('launch-view-run').style.display     = '';
        document.getElementById('launch-view-analysis').style.display = 'none';
        document.getElementById('launch-view-advanced').style.display = 'none';
        document.getElementById('launch-view-docker').style.display = 'none';
        document.getElementById('launch-tab-run').classList.add('active');
        document.getElementById('launch-tab-analysis').classList.remove('active');
        document.getElementById('launch-tab-advanced').classList.remove('active');
        document.getElementById('launch-tab-docker').classList.remove('active');
        document.getElementById('launch-result').style.display = 'none';
    });
    document.getElementById('launch-tab-analysis')?.addEventListener('click', () => {
        activeLaunchSubTab = 'analysis';
        document.getElementById('launch-view-run').style.display     = 'none';
        document.getElementById('launch-view-analysis').style.display = '';
        document.getElementById('launch-view-advanced').style.display = 'none';
        document.getElementById('launch-view-docker').style.display = 'none';
        document.getElementById('launch-tab-run').classList.remove('active');
        document.getElementById('launch-tab-analysis').classList.add('active');
        document.getElementById('launch-tab-advanced').classList.remove('active');
        document.getElementById('launch-tab-docker').classList.remove('active');
        document.getElementById('launch-result').style.display = 'none';
    });
    document.getElementById('launch-tab-advanced')?.addEventListener('click', () => {
        activeLaunchSubTab = 'advanced';
        document.getElementById('launch-view-run').style.display     = 'none';
        document.getElementById('launch-view-analysis').style.display = 'none';
        document.getElementById('launch-view-advanced').style.display = '';
        document.getElementById('launch-view-docker').style.display = 'none';
        document.getElementById('launch-tab-run').classList.remove('active');
        document.getElementById('launch-tab-analysis').classList.remove('active');
        document.getElementById('launch-tab-advanced').classList.add('active');
        document.getElementById('launch-tab-docker').classList.remove('active');
        document.getElementById('launch-result').style.display = 'none';
    });
    document.getElementById('launch-tab-docker')?.addEventListener('click', () => {
        activeLaunchSubTab = 'docker';
        document.getElementById('launch-view-run').style.display     = 'none';
        document.getElementById('launch-view-analysis').style.display = 'none';
        document.getElementById('launch-view-advanced').style.display = 'none';
        document.getElementById('launch-view-docker').style.display = '';
        document.getElementById('launch-tab-run').classList.remove('active');
        document.getElementById('launch-tab-analysis').classList.remove('active');
        document.getElementById('launch-tab-advanced').classList.remove('active');
        document.getElementById('launch-tab-docker').classList.add('active');
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

    // Convert a JSON object to a CLI string: {run:"X", sample_filter:["S1","S2"]} -> "--run=X --sample_filter=S1,S2"
    function jsonCommandToCli(obj) {
        return Object.entries(obj).map(([k, v]) => {
            if (Array.isArray(v))           return `--${k}=${v.join(',')}`;
            if (v === true)                 return `--${k}`;
            if (v === false || v === null)  return '';
            return `--${k}=${v}`;
        }).filter(Boolean).join(' ');
    }

    function applyAnalysisModuleDefaults() {
        const minput = document.getElementById('analysis-module');
        if (!minput) return;
        const meta = (window._modulesMeta || {})[minput.value.toUpperCase()];
        if (!meta) return;
        const defaults = meta.defaults || {};
        const qsel = document.getElementById('analysis-queue');
        if (qsel && defaults.queue) qsel.value = defaults.queue;
        const thr = document.getElementById('analysis-threads');
        if (thr && defaults.threads) thr.value = defaults.threads;
        const mem = document.getElementById('analysis-memory');
        if (mem && defaults.memory) mem.value = defaults.memory;
        const pri = document.getElementById('analysis-prioritize');
        if (pri) pri.checked = !!defaults.prioritize;
    }

    // Switch command input panels when user changes format
    document.querySelectorAll('input[name="analysis-cmd-mode"]').forEach(radio => {
        radio.addEventListener('change', () => {
            const mode = document.querySelector('input[name="analysis-cmd-mode"]:checked')?.value;
            document.getElementById('analysis-cmd-text').style.display = mode === 'text' ? '' : 'none';
            document.getElementById('analysis-cmd-json').style.display = mode === 'json' ? '' : 'none';
            document.getElementById('analysis-cmd-file').style.display = mode === 'file' ? '' : 'none';
        });
    });

    document.getElementById('analysis-form-analysis')?.addEventListener('submit', async e => {
        e.preventDefault();
        const msel = document.getElementById('analysis-module');
        const moduleName = msel?.value || '';

        const mode = document.querySelector('input[name="analysis-cmd-mode"]:checked')?.value || 'text';
        let command = '';
        if (mode === 'text') {
            command = document.getElementById('analysis-command-text')?.value.trim() || '';
        } else if (mode === 'json') {
            const raw = document.getElementById('analysis-command-json')?.value.trim() || '';
            try { command = jsonCommandToCli(JSON.parse(raw)); }
            catch (_) { showLaunchResult('Invalid JSON command', false); return; }
        } else if (mode === 'file') {
            const file = document.getElementById('analysis-command-file')?.files[0];
            if (!file) { showLaunchResult('No file selected', false); return; }
            try {
                const text = await file.text();
                command = jsonCommandToCli(JSON.parse(text));
            } catch (_) { showLaunchResult('Invalid JSON file', false); return; }
        }
        if (!command) { showLaunchResult('Command is required', false); return; }

        const analysis_name = document.getElementById('analysis-name-input')?.value.trim() || undefined;
        const queue = document.getElementById('analysis-queue')?.value || undefined;
        const threadsRaw = document.getElementById('analysis-threads')?.value.trim();
        const threads = threadsRaw ? parseInt(threadsRaw, 10) : undefined;
        const memory = document.getElementById('analysis-memory')?.value.trim() || undefined;
        const prioritize = document.getElementById('analysis-prioritize')?.checked || undefined;

        const payload = { module: moduleName, command };
        if (analysis_name) payload.analysis_name = analysis_name;
        if (queue)         payload.queue = queue;
        if (threads)       payload.threads = threads;
        if (memory)        payload.memory = memory;
        if (prioritize)    payload.prioritize = prioritize;

        const resp = await fetch('/analysis', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${token}` },
            body: JSON.stringify(payload),
        });
        const text = await resp.text();
        showLaunchResult(text, resp.ok);
        if (resp.ok) setTimeout(fetchAllTasks, 500);
    });

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

    document.getElementById('analysis-form-docker')?.addEventListener('submit', async e => {
        e.preventDefault();
        const analysis_name = document.getElementById('docker-analysis-name')?.value.trim();
        if (!analysis_name) {
            showLaunchResult('Analysis name is required', false);
            return;
        }
        const image = document.getElementById('docker-image')?.value.trim();
        if (!analysis_name || !image) {
            showLaunchResult('Docker image is required', false);
            return;
        }
        const command = document.getElementById('docker-command')?.value.trim();
        if (!command) {
            showLaunchResult('Docker command is required', false);
            return;
        }
        const docker_extra_params = document.getElementById('docker-extra-params')?.value.trim();
        const use_stark_container_mount = document.getElementById('docker-use-stark-container-mount')?.checked;
        const queue = document.getElementById('docker-queue')?.value.trim();
        const threads = Number.parseInt(document.getElementById('docker-threads')?.value.trim());
        const memory = document.getElementById('docker-memory')?.value.trim();
        const prioritize = document.getElementById('docker-prioritize')?.checked;
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
            archivesStatusDd = buildDropdown('archives-filter-status-container', statuses, null, renderArchivesTable, 'All states');
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
        const statusCls = { finished: 'state-finished-success', failed: 'state-finished-failed', unknown: 'state-unknown' };
        const table = document.createElement('table');
        table.className = 'cluster-tasks-table';
        table.innerHTML = `<thead><tr>
            <th data-col="node">Node <span class="sort-icon"></span></th>
            <th data-col="queue">Queue <span class="sort-icon"></span></th>
            <th data-col="threads">Slots <span class="sort-icon"></span></th>
            <th data-col="status">State <span class="sort-icon"></span></th>
            <th data-col="mtime">Time <span class="sort-icon"></span></th>
            <th data-col="run_name">Analysis Name <span class="sort-icon"></span></th>
            <th>Actions</th>
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
            let date_title = a.mtime != null ? `Launch:\t${formatDate(a.mtime) || '-'}` : 'Date unknown';
            date_title += a.end_date != null ? `\nEnd:\t\t${formatDate(a.end_date) || '-'}` : '';
            date_title += a.exec_time != null ? `\nTime:\t${formatTime(a.exec_time) || '-'}` : '';
            const analysis_name = a.run_name ? a.run_name : '-';
            const analysis_name_label = a.run_name ? (a.run_name.length > 50 ? a.run_name.slice(0, 47) + '...' : a.run_name) : '-';
            const analysis_id = a.analysis_id_name;
            const analysis_title = a.analysis_id_name ? `Analysis Name:\t${analysis_name}\nAnalysis ID:\t${analysis_id}` : '';
            const node_title = a.node ? `Node:\t${a.node}\nURL:\t\t${a.node_url || '-'}` : 'Node unknown';
            const node_label = a.node ? (a.node.length > 12 ? a.node.slice(0, 12) + '...' : a.node) : '-';
            tr.innerHTML = `
                <td title="${node_title}">${node_label}</td>
                <td>${a.queue || '-'}</td>
                <td>${a.threads ?? '-'}</td>
                <td class="${statusCls[statusLabel] || 'state-unknown'}">${statusLabel}</td>
                <td title="${date_title}">${formatDate(a.mtime) || '-'}</td>
                <td title="${analysis_title}" >${analysis_name_label}</td>
            `;
            const id = a.analysis_id_name;
            const nodeUrl = a.node_url || '';

            const actionTd = document.createElement('td');
            actionTd.className = 'action-buttons';

            // L — Log, A — Analysis JSON (read-only, no admin required)
            // R — Relaunch, X — Delete (admin only)
            const btnDefs = [['L', 'log'], ['A', 'json'], ['R', 'relaunch'], ['X', 'delete']];
            btnDefs.forEach(([label, act]) => {
                const isDanger = ['relaunch', 'delete'].includes(act);
                if (isDanger && !isAdmin) return;
                const btn = document.createElement('button');
                btn.textContent = label;
                btn.title = `${capitalize(act)} '${a.run_name || id}'`;
                if (isDanger) btn.classList.add('btn-danger');

                if (act === 'log' || act === 'json') {
                    btn.addEventListener('click', () => {
                        const key = `archive::${id}::${act}`;
                        const params = new URLSearchParams({ node_url: nodeUrl });
                        toggleDetail(archivesOpenDetails, tr, `/cluster/proxy/archive/${encodeURIComponent(id)}/${act}?${params}`, key);
                    });
                } else {
                    btn.addEventListener('click', async () => {
                        const actionLabel = act === 'delete'
                            ? `Delete all files for "${a.run_name || id}"?\n\nThis will permanently remove the .json, .info and .output files.`
                            : `Relaunch '${a.run_name || id}'?`;
                        if (!confirm(actionLabel)) return;
                        btn.disabled = true; btn.textContent = '...';
                        let ok = false;
                        try {
                            const params = new URLSearchParams({ node_url: nodeUrl });
                            let resp;
                            if (act === 'relaunch') {
                                resp = await fetch(`/cluster/proxy/archive/${encodeURIComponent(id)}/relaunch?${params}`, {
                                    method: 'POST',
                                    headers: { Authorization: `Bearer ${token}` },
                                });
                            } else {
                                resp = await fetch(`/cluster/proxy/archive/${encodeURIComponent(id)}?${params}`, {
                                    method: 'DELETE',
                                    headers: { Authorization: `Bearer ${token}` },
                                });
                            }
                            ok = resp.ok;
                            if (ok && act === 'delete') {
                                allArchives = allArchives.filter(x => x.analysis_id_name !== id);
                            }
                        } catch (_) {}
                        btn.textContent = ok ? '\u2713 Done' : '\u2717 Failed';
                        btn.classList.toggle('btn-success', ok);
                        btn.classList.toggle('btn-error',   !ok);
                        setTimeout(() => {
                            btn.textContent = label; btn.disabled = false;
                            btn.classList.remove('btn-success', 'btn-error');
                            if (act === 'delete' && ok) renderArchivesTable();
                            else if (act === 'relaunch' && ok) fetchAllTasks();
                        }, 1500);
                    });
                }
                actionTd.appendChild(btn);
            });

            tr.appendChild(actionTd);
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
        el.innerHTML = '<p class="cluster-loading">Loading...</p>';
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

    // ══════════════════════════════════════════════════════════════════════════
    // STATISTICS TAB
    // ══════════════════════════════════════════════════════════════════════════

    let statsBarChart     = null;
    let statsDonutChart   = null;
    let statsQueueDd      = null;
    let statsNodeDd       = null;
    let activeStatsSubTab = 'distribution';
    let currentStats      = null;

    const STATS_WARN_THRESHOLD = 90;

    const STATS_PRESETS = {
        year:  [{v:'3y',l:'Last 3 years'},{v:'5y',l:'Last 5 years'},{v:'10y',l:'Last 10 years'},{v:'all',l:'All'},{v:'custom',l:'Custom...'}],
        month: [{v:'3m',l:'Last 3 months'},{v:'6m',l:'Last 6 months'},{v:'12m',l:'Last 12 months'},{v:'24m',l:'Last 24 months'},{v:'all',l:'All'},{v:'custom',l:'Custom...'}],
        day:   [{v:'7d',l:'Last 7 days'},{v:'14d',l:'Last 14 days'},{v:'30d',l:'Last 30 days'},{v:'90d',l:'Last 90 days'},{v:'all',l:'All'},{v:'custom',l:'Custom...'}],
    };

    function updateStatsPresets(keepValue = false) {
        const granularity = document.getElementById('stats-granularity')?.value || 'month';
        const sel = document.getElementById('stats-daterange-preset');
        if (!sel) return;
        const presets = STATS_PRESETS[granularity] || STATS_PRESETS.month;
        const current = sel.value;
        sel.innerHTML = '';
        presets.forEach(p => {
            const opt = document.createElement('option');
            opt.value = p.v; opt.textContent = p.l;
            sel.appendChild(opt);
        });
        if (keepValue && presets.some(p => p.v === current)) sel.value = current;
        const customGroup = document.getElementById('stats-custom-range-group');
        if (customGroup) customGroup.style.display = sel.value === 'custom' ? 'flex' : 'none';
    }

    function getStatsDateRange() {
        const preset = document.getElementById('stats-daterange-preset')?.value || 'all';
        if (preset === 'all') return { from: null, to: null };
        if (preset === 'custom') {
            const fromVal = document.getElementById('stats-date-from')?.value;
            const toVal   = document.getElementById('stats-date-to')?.value;
            return {
                from: fromVal ? new Date(fromVal).getTime() / 1000         : null,
                to:   toVal   ? new Date(toVal).getTime()   / 1000 + 86399 : null,
            };
        }
        const m = preset.match(/^(\d+)([dmy])$/);
        if (!m) return { from: null, to: null };
        const n = parseInt(m[1]), unit = m[2];
        const d = new Date();
        if      (unit === 'd') d.setDate(d.getDate() - n);
        else if (unit === 'm') d.setMonth(d.getMonth() - n);
        else                   d.setFullYear(d.getFullYear() - n);
        return { from: d.getTime() / 1000, to: null };
    }

    function showStatsWarning(labelCount) {
        const el = document.getElementById('stats-warning');
        if (!el) return;
        if (labelCount > STATS_WARN_THRESHOLD) {
            el.style.display = '';
            el.textContent = `\u26a0\ufe0f ${labelCount} periods to display \u2014 the chart may be hard to read. Consider switching to a larger granularity (Month or Year) or narrowing the date range.`;
        } else {
            el.style.display = 'none';
        }
    }

    function buildStatsQueueFilter(archives) {
        const queues = [...new Set(archives.map(a => a.queue || ''))].filter(Boolean).sort();
        const group  = document.getElementById('stats-filter-queues-group');
        if (queues.length > 1) {
            if (group) group.style.display = 'flex';
            if (!statsQueueDd)
                statsQueueDd = buildDropdown('stats-filter-queues-container', queues, null, renderStats, 'All queues');
        }
        const nodes     = [...new Set(archives.map(a => a.node || ''))].filter(Boolean).sort();
        const nodeGroup = document.getElementById('stats-filter-nodes-group');
        if (nodes.length > 1) {
            if (nodeGroup) nodeGroup.style.display = 'flex';
            if (!statsNodeDd)
                statsNodeDd = buildDropdown('stats-filter-nodes-container', nodes, null, renderStats, 'All nodes');
        }
    }

    function getPeriodKey(ts, granularity) {
        if (ts == null) return null;
        const d   = new Date(ts * 1000);
        const y   = d.getFullYear();
        const mo  = String(d.getMonth() + 1).padStart(2, '0');
        const day = String(d.getDate()).padStart(2, '0');
        if (granularity === 'year')  return `${y}`;
        if (granularity === 'month') return `${y}-${mo}`;
        return `${y}-${mo}-${day}`;
    }

    function computeStats(archives, granularity, metric, dateRange) {
        const selectedQueues = statsQueueDd ? statsQueueDd.getSelected() : null;
        const selectedNodes  = statsNodeDd  ? statsNodeDd.getSelected()  : null;
        let filtered = selectedQueues
            ? archives.filter(a => selectedQueues.includes(a.queue || ''))
            : [...archives];
        if (selectedNodes) filtered = filtered.filter(a => selectedNodes.includes(a.node || ''));
        if (dateRange.from != null) filtered = filtered.filter(a => a.mtime >= dateRange.from);
        if (dateRange.to   != null) filtered = filtered.filter(a => a.mtime <= dateRange.to);

        const periodMap = new Map();
        for (const a of filtered) {
            const key = getPeriodKey(a.mtime, granularity);
            if (key === null) continue;
            if (!periodMap.has(key)) periodMap.set(key, { finished: 0, failed: 0, unknown: 0, total: 0 });
            const entry  = periodMap.get(key);
            const weight = metric === 'slots' ? (a.threads ?? 1) : 1;
            const status = a.status || 'unknown';
            if (status === 'finished')    entry.finished += weight;
            else if (status === 'failed') entry.failed   += weight;
            else                          entry.unknown  += weight;
            entry.total += weight;
        }

        const labels       = [...periodMap.keys()].sort();
        const finishedData = labels.map(l => periodMap.get(l).finished);
        const failedData   = labels.map(l => periodMap.get(l).failed);
        const unknownData  = labels.map(l => periodMap.get(l).unknown);

        let totalFinished = 0, totalFailed = 0, totalUnknown = 0;
        for (const e of periodMap.values()) {
            totalFinished += e.finished;
            totalFailed   += e.failed;
            totalUnknown  += e.unknown;
        }

        const tableRows = labels.map(l => {
            const e   = periodMap.get(l);
            const rate = e.total > 0 ? Math.round(e.finished / e.total * 100) : 0;
            return { period: l, finished: e.finished, failed: e.failed, unknown: e.unknown, total: e.total, rate };
        });

        return { labels, finishedData, failedData, unknownData, totalFinished, totalFailed, totalUnknown, tableRows };
    }

    function renderStatsBarChart(stats, metric) {
        const metricLabel = metric === 'slots' ? 'Slots' : 'Tasks';
        const barCanvas = document.getElementById('stats-bar-chart');
        const noData    = document.getElementById('stats-bar-no-data');
        if (!barCanvas || typeof Chart === 'undefined') return;
        if (statsBarChart) { statsBarChart.destroy(); statsBarChart = null; }
        if (!stats.labels.length) {
            barCanvas.style.display = 'none';
            if (noData) { noData.style.display = ''; noData.textContent = 'No data in selected range.'; }
            return;
        }
        barCanvas.style.display = '';
        if (noData) noData.style.display = 'none';
        statsBarChart = new Chart(barCanvas, {
            type: 'bar',
            data: {
                labels: stats.labels,
                datasets: [
                    { label: 'Finished', data: stats.finishedData, backgroundColor: 'rgba(40,167,69,0.82)',   stack: 's' },
                    { label: 'Failed',   data: stats.failedData,   backgroundColor: 'rgba(220,53,69,0.82)',   stack: 's' },
                    { label: 'Unknown',  data: stats.unknownData,  backgroundColor: 'rgba(108,117,125,0.45)', stack: 's' },
                ],
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: { position: 'top' },
                    tooltip: { mode: 'index', intersect: false },
                },
                scales: {
                    x: { stacked: true, ticks: { maxRotation: 45, autoSkip: true, maxTicksLimit: 60 } },
                    y: { stacked: true, beginAtZero: true, title: { display: true, text: metricLabel } },
                },
            },
        });
    }

    function renderStatsDonutChart(stats) {
        const donutCanvas = document.getElementById('stats-donut-chart');
        const noData      = document.getElementById('stats-donut-no-data');
        if (!donutCanvas || typeof Chart === 'undefined') return;
        if (statsDonutChart) { statsDonutChart.destroy(); statsDonutChart = null; }
        const grand = stats.totalFinished + stats.totalFailed + stats.totalUnknown;
        if (!grand) {
            donutCanvas.style.display = 'none';
            if (noData) { noData.style.display = ''; noData.textContent = 'No data in selected range.'; }
            return;
        }
        donutCanvas.style.display = '';
        if (noData) noData.style.display = 'none';
        const rate      = Math.round(stats.totalFinished / grand * 100);
        const rateColor = rate >= 80 ? '#28a745' : rate >= 50 ? '#e07b00' : '#dc3545';
        statsDonutChart = new Chart(donutCanvas, {
            type: 'doughnut',
            data: {
                labels: ['Finished', 'Failed', 'Unknown'],
                datasets: [{
                    data: [stats.totalFinished, stats.totalFailed, stats.totalUnknown],
                    backgroundColor: ['rgba(40,167,69,0.82)', 'rgba(220,53,69,0.82)', 'rgba(108,117,125,0.45)'],
                    borderColor:     ['#28a745',               '#dc3545',               '#6c757d'],
                    borderWidth: 1,
                }],
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                cutout: '62%',
                plugins: {
                    legend: { position: 'bottom' },
                    tooltip: {
                        callbacks: {
                            label: ctx => {
                                const v   = ctx.raw;
                                const pct = grand > 0 ? Math.round(v / grand * 100) : 0;
                                return ` ${ctx.label}: ${v} (${pct}%)`;
                            },
                        },
                    },
                },
            },
            plugins: [{
                id: 'statsCenterText',
                afterDraw(chart) {
                    const { ctx, chartArea: { top, bottom, left, right } } = chart;
                    const cx = (left + right) / 2;
                    const cy = (top  + bottom) / 2;
                    ctx.save();
                    ctx.font         = 'bold 1.8em sans-serif';
                    ctx.fillStyle    = rateColor;
                    ctx.textAlign    = 'center';
                    ctx.textBaseline = 'middle';
                    ctx.fillText(`${rate}%`, cx, cy - 10);
                    ctx.font      = '0.78em sans-serif';
                    ctx.fillStyle = '#888';
                    ctx.fillText('success', cx, cy + 14);
                    ctx.restore();
                },
            }],
        });
    }

    function renderStatsTable(stats, metric) {
        const el = document.getElementById('stats-table-content');
        if (!el) return;
        const metricLabel = metric === 'slots' ? 'Slots' : 'Tasks';
        if (!stats.tableRows.length) {
            el.innerHTML = '<p style="color:#888;font-style:italic">No data in selected range.</p>'; return;
        }
        const table = document.createElement('table');
        table.className = 'cluster-tasks-table stats-summary-table';
        table.innerHTML = `<thead><tr>
            <th>Period</th>
            <th>Finished</th>
            <th>Failed</th>
            <th>Unknown</th>
            <th>Total ${metricLabel}</th>
            <th>Success Rate</th>
        </tr></thead>`;
        const tbody = document.createElement('tbody');
        stats.tableRows.forEach(r => {
            const cls = r.rate >= 80 ? 'stats-rate-good' : r.rate >= 50 ? 'stats-rate-warn' : 'stats-rate-bad';
            const tr  = document.createElement('tr');
            tr.innerHTML = `
                <td>${r.period}</td>
                <td class="state-finished-success">${r.finished}</td>
                <td class="state-finished-failed">${r.failed}</td>
                <td style="color:#6c757d">${r.unknown}</td>
                <td>${r.total}</td>
                <td><span class="stats-rate-badge ${cls}">${r.rate}%</span></td>`;
            tbody.appendChild(tr);
        });
        const tot = stats.tableRows.reduce(
            (acc, r) => ({ finished: acc.finished + r.finished, failed: acc.failed + r.failed, unknown: acc.unknown + r.unknown, total: acc.total + r.total }),
            { finished: 0, failed: 0, unknown: 0, total: 0 }
        );
        const totalRate = tot.total > 0 ? Math.round(tot.finished / tot.total * 100) : 0;
        const totalCls  = totalRate >= 80 ? 'stats-rate-good' : totalRate >= 50 ? 'stats-rate-warn' : 'stats-rate-bad';
        const totTr = document.createElement('tr');
        totTr.className = 'stats-total-row';
        totTr.innerHTML = `
            <td><strong>TOTAL</strong></td>
            <td class="state-finished-success"><strong>${tot.finished}</strong></td>
            <td class="state-finished-failed"><strong>${tot.failed}</strong></td>
            <td style="color:#6c757d"><strong>${tot.unknown}</strong></td>
            <td><strong>${tot.total}</strong></td>
            <td><span class="stats-rate-badge ${totalCls}">${totalRate}%</span></td>`;
        tbody.appendChild(totTr);
        table.appendChild(tbody);
        el.innerHTML = ''; el.appendChild(table);
    }

    function switchStatsSubTab(name) {
        activeStatsSubTab = name;
        ['distribution', 'summary', 'success'].forEach(n => {
            const view = document.getElementById(`stats-view-${n}`);
            const btn  = document.getElementById(`stats-tab-${n}`);
            if (view) view.style.display = n === name ? '' : 'none';
            if (btn)  btn.classList.toggle('active', n === name);
        });
        renderActiveStatsSubTab();
    }

    function renderActiveStatsSubTab() {
        if (!currentStats) return;
        const metric = document.getElementById('stats-metric')?.value || 'tasks';
        if      (activeStatsSubTab === 'distribution') renderStatsBarChart(currentStats, metric);
        else if (activeStatsSubTab === 'summary')      renderStatsTable(currentStats, metric);
        else if (activeStatsSubTab === 'success')      renderStatsDonutChart(currentStats);
    }

    function renderStats() {
        if (!allArchives.length) return;
        const granularity = document.getElementById('stats-granularity')?.value || 'month';
        const metric      = document.getElementById('stats-metric')?.value      || 'tasks';
        const dateRange   = getStatsDateRange();
        currentStats      = computeStats(allArchives, granularity, metric, dateRange);
        showStatsWarning(currentStats.labels.length);
        renderActiveStatsSubTab();
    }

    async function fetchStats() {
        const el = document.getElementById('stats-table-content');
        if (!allArchives.length) {
            if (el) el.innerHTML = '<p class="cluster-loading">Loading...</p>';
            try {
                const resp = await fetch('/cluster/archives', { headers: { Authorization: `Bearer ${token}` } });
                if (!resp.ok) { if (el) el.innerHTML = `<p class="error">Error ${resp.status}</p>`; return; }
                const raw = await resp.json();
                const seen = new Set();
                allArchives = raw.filter(a => {
                    const key = a.analysis_id_name;
                    if (seen.has(key)) return false;
                    seen.add(key);
                    return true;
                });
            } catch (e) {
                if (el) el.innerHTML = `<p class="error">Could not load archives: ${e}</p>`;
                return;
            }
        }
        buildStatsQueueFilter(allArchives);
        updateStatsPresets(true);
        renderStats();
    }

    document.getElementById('stats-refresh-btn')?.addEventListener('click', () => {
        allArchives  = [];
        statsQueueDd = null;
        statsNodeDd  = null;
        currentStats = null;
        const grp = document.getElementById('stats-filter-queues-group');
        if (grp) { grp.style.display = 'none'; const c = document.getElementById('stats-filter-queues-container'); if (c) c.innerHTML = ''; }
        const ngrp = document.getElementById('stats-filter-nodes-group');
        if (ngrp) { ngrp.style.display = 'none'; const c = document.getElementById('stats-filter-nodes-container'); if (c) c.innerHTML = ''; }
        fetchStats();
    });
    document.getElementById('stats-granularity')?.addEventListener('change', () => { updateStatsPresets(false); renderStats(); });
    document.getElementById('stats-metric')?.addEventListener('change', renderStats);
    document.getElementById('stats-daterange-preset')?.addEventListener('change', () => {
        const val = document.getElementById('stats-daterange-preset')?.value;
        const customGroup = document.getElementById('stats-custom-range-group');
        if (customGroup) customGroup.style.display = val === 'custom' ? 'flex' : 'none';
        if (val !== 'custom') renderStats();
    });
    document.getElementById('stats-date-from')?.addEventListener('change', renderStats);
    document.getElementById('stats-date-to')?.addEventListener('change', renderStats);
    ['distribution', 'summary', 'success'].forEach(name => {
        document.getElementById(`stats-tab-${name}`)?.addEventListener('click', () => switchStatsSubTab(name));
    });
    updateStatsPresets();

    // ── Auto-refresh ──────────────────────────────────────────────────────────
    const REFRESH_MS = typeof REFRESH_INTERVAL_MS !== 'undefined' ? REFRESH_INTERVAL_MS : 10000;
    setInterval(() => {
        if (activeTab === 'tasks')   fetchAllTasks();
        if (activeTab === 'cluster') fetchClusterSummary();
    }, REFRESH_MS);

});
