import "core-js/stable";
import 'bootstrap';
import './pdb-redo-style.scss';

class Tortoize {
	constructor(form) {
		this.form = document.forms['tortoize-form'];
		this.resultDiv = document.getElementById('tortoize-result');
		this.table = document.getElementById('tortoize-table');
		this.alert = document.getElementById('tortoize-alert');
		this.form.addEventListener('submit', (evt) => this.handleSubmit(evt));

		this.alert.classList.add('invisible');
	}

	async handleSubmit(event) {
		if (event)
			event.preventDefault();

		try {
			const fd = new FormData(this.form);
			this.resultDiv.style.display = 'none';
			this.alert.classList.add('invisible');

			const r = await fetch(this.form.action, {
				credentials: 'include',
				method: 'post',
				body: fd
			})

			const data = await r.json()
			if (typeof data.error === "string")
				throw data.error;

			if (!r.ok)
				throw `Fetching result failed, ${r.statusText}`;

			if (data.model) {
				this.process(data.model);

				const downloadBtn = document.getElementById("download-btn");
				downloadBtn.addEventListener("click", (e) => {
					//convert JSON data to URI data scheme.
					const url = "data:application/json;charset=utf-8," + encodeURIComponent(JSON.stringify(data));

					const link = document.createElement("a");
					link.href = url;
					link.download = "tortoize-result.json";
					document.body.appendChild(link);
					link.click();
					document.body.removeChild(link);
				});
			}
			else
				throw 'Reply does not contain data';
		}

		catch (err) {
			console.log(err);

			this.alert.textContent = `Could not calculate rama z-score: ${err}`;
			this.alert.classList.remove('invisible');
		}
	}

	process(model) {
		const tbody = this.table.querySelector('tbody');
		[...tbody.querySelectorAll('tr')].forEach(row => tbody.removeChild(row));

		for (const [id, data] of Object.entries(model)) {
			const row = document.createElement("tr");
			const tdid = document.createElement("td");
			tdid.textContent = id;
			row.appendChild(tdid);

			for (const f of ['ramachandran', 'torsion']) {
				const td = document.createElement('td');
				td.textContent = (+data[`${f}-z`]).toFixed(2) + " ± " + (+data[`${f}-jackknife-sd`]).toFixed(2);
				row.appendChild(td);
			}

			tbody.appendChild(row);
		}


		this.resultDiv.style.display = 'unset';
	}
}


window.addEventListener('load', () => {
	new Tortoize();
});